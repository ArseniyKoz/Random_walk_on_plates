#include "wop/config/config.hpp"

#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace wop::config {

namespace {

struct Node {
    enum class Type {
        Scalar,
        Map,
        Sequence,
    };

    using MapType = std::map<std::string, Node>;
    using SequenceType = std::vector<Node>;

    Type type = Type::Scalar;
    std::string scalar;
    MapType map;
    SequenceType sequence;

    static Node make_scalar(std::string value) {
        Node node;
        node.type = Type::Scalar;
        node.scalar = std::move(value);
        return node;
    }

    static Node make_map(MapType value) {
        Node node;
        node.type = Type::Map;
        node.map = std::move(value);
        return node;
    }

    static Node make_sequence(SequenceType value) {
        Node node;
        node.type = Type::Sequence;
        node.sequence = std::move(value);
        return node;
    }
};

struct Line {
    int indent = 0;
    std::string text;
};

std::string trim(const std::string& text) {
    std::size_t begin = 0;
    while (begin < text.size() && std::isspace(static_cast<unsigned char>(text[begin])) != 0) {
        ++begin;
    }
    std::size_t end = text.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
        --end;
    }
    return text.substr(begin, end - begin);
}

std::string strip_comments(std::string text) {
    bool in_single = false;
    bool in_double = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (ch == '\'' && !in_double) {
            in_single = !in_single;
            continue;
        }
        if (ch == '"' && !in_single) {
            in_double = !in_double;
            continue;
        }
        if (ch == '#' && !in_single && !in_double) {
            text.erase(i);
            break;
        }
    }
    return text;
}

std::vector<Line> preprocess_lines(const std::string& raw_text) {
    std::vector<Line> lines;
    std::istringstream input(raw_text);
    std::string raw_line;
    while (std::getline(input, raw_line)) {
        if (!raw_line.empty() && raw_line.back() == '\r') {
            raw_line.pop_back();
        }
        raw_line = strip_comments(raw_line);
        if (trim(raw_line).empty()) {
            continue;
        }
        if (raw_line.find('\t') != std::string::npos) {
            throw std::invalid_argument("Tabs are not supported in YAML config.");
        }

        int indent = 0;
        while (indent < static_cast<int>(raw_line.size()) && raw_line[static_cast<std::size_t>(indent)] == ' ') {
            ++indent;
        }
        if ((indent % 2) != 0) {
            throw std::invalid_argument("YAML indentation must use multiples of two spaces.");
        }
        lines.push_back(Line{indent, raw_line.substr(static_cast<std::size_t>(indent))});
    }
    return lines;
}

std::pair<std::string, std::string> split_key_value(const std::string& text) {
    const std::size_t pos = text.find(':');
    if (pos == std::string::npos) {
        throw std::invalid_argument("Expected key: value entry in YAML config.");
    }
    return {trim(text.substr(0, pos)), trim(text.substr(pos + 1))};
}

Node parse_block(const std::vector<Line>& lines, std::size_t& index, int indent);

void parse_map_entry(
    const std::vector<Line>& lines,
    std::size_t& index,
    int indent,
    std::map<std::string, Node>& out_map,
    const std::string& text) {
    const auto [key, value] = split_key_value(text);
    if (key.empty()) {
        throw std::invalid_argument("YAML key must not be empty.");
    }

    if (!value.empty()) {
        out_map[key] = Node::make_scalar(value);
        ++index;
        return;
    }

    ++index;
    if (index >= lines.size() || lines[index].indent <= indent) {
        throw std::invalid_argument("Expected nested YAML block after key: " + key);
    }
    out_map[key] = parse_block(lines, index, indent + 2);
}

Node parse_map(const std::vector<Line>& lines, std::size_t& index, int indent) {
    std::map<std::string, Node> result;
    while (index < lines.size()) {
        if (lines[index].indent < indent) {
            break;
        }
        if (lines[index].indent > indent) {
            throw std::invalid_argument("Unexpected indentation in YAML map.");
        }
        if (lines[index].text.rfind("- ", 0) == 0) {
            break;
        }
        parse_map_entry(lines, index, indent, result, lines[index].text);
    }
    return Node::make_map(std::move(result));
}

Node parse_sequence(const std::vector<Line>& lines, std::size_t& index, int indent) {
    std::vector<Node> result;
    while (index < lines.size() && lines[index].indent == indent && lines[index].text.rfind("- ", 0) == 0) {
        const std::string rest = trim(lines[index].text.substr(2));
        if (rest.empty()) {
            ++index;
            if (index >= lines.size() || lines[index].indent <= indent) {
                throw std::invalid_argument("Expected nested YAML block after sequence item.");
            }
            result.push_back(parse_block(lines, index, indent + 2));
            continue;
        }

        if (rest.find(':') != std::string::npos) {
            std::map<std::string, Node> item_map;
            parse_map_entry(lines, index, indent, item_map, rest);
            while (index < lines.size() && lines[index].indent == indent + 2 && lines[index].text.rfind("- ", 0) != 0) {
                parse_map_entry(lines, index, indent + 2, item_map, lines[index].text);
            }
            result.push_back(Node::make_map(std::move(item_map)));
            continue;
        }

        result.push_back(Node::make_scalar(rest));
        ++index;
    }
    return Node::make_sequence(std::move(result));
}

Node parse_block(const std::vector<Line>& lines, std::size_t& index, int indent) {
    if (index >= lines.size() || lines[index].indent != indent) {
        throw std::invalid_argument("YAML block indentation mismatch.");
    }
    if (lines[index].text.rfind("- ", 0) == 0) {
        return parse_sequence(lines, index, indent);
    }
    return parse_map(lines, index, indent);
}

Node parse_yaml_subset(const std::string& raw_text) {
    const auto lines = preprocess_lines(raw_text);
    if (lines.empty()) {
        throw std::invalid_argument("Config file is empty.");
    }
    std::size_t index = 0;
    Node root = parse_block(lines, index, 0);
    if (index != lines.size()) {
        throw std::invalid_argument("Could not parse entire YAML config.");
    }
    return root;
}

const Node& require_map(const Node& node, const std::string& name) {
    if (node.type != Node::Type::Map) {
        throw std::invalid_argument(name + " must be a YAML mapping.");
    }
    return node;
}

const Node& require_sequence(const Node& node, const std::string& name) {
    if (node.type != Node::Type::Sequence) {
        throw std::invalid_argument(name + " must be a YAML sequence.");
    }
    return node;
}

const std::string& require_scalar(const Node& node, const std::string& name) {
    if (node.type != Node::Type::Scalar) {
        throw std::invalid_argument(name + " must be a YAML scalar.");
    }
    return node.scalar;
}

const Node* find_key(const Node& map_node, const std::string& key) {
    if (map_node.type != Node::Type::Map) {
        return nullptr;
    }
    const auto it = map_node.map.find(key);
    return (it == map_node.map.end()) ? nullptr : &it->second;
}

const Node& require_key(const Node& map_node, const std::string& key, const std::string& parent_name) {
    const Node* node = find_key(map_node, key);
    if (node == nullptr) {
        throw std::invalid_argument("Missing required key " + parent_name + "." + key + ".");
    }
    return *node;
}

std::string unquote(std::string text) {
    text = trim(text);
    if (text.size() >= 2) {
        const char first = text.front();
        const char last = text.back();
        if ((first == '"' && last == '"') || (first == '\'' && last == '\'')) {
            return text.substr(1, text.size() - 2);
        }
    }
    return text;
}

double parse_double_scalar(const std::string& text, const std::string& name) {
    const std::string clean = unquote(text);
    std::size_t consumed = 0;
    const double value = std::stod(clean, &consumed);
    if (consumed != clean.size()) {
        throw std::invalid_argument("Invalid floating-point value for " + name + ".");
    }
    return value;
}

int parse_int_scalar(const std::string& text, const std::string& name) {
    const std::string clean = unquote(text);
    std::size_t consumed = 0;
    const int value = std::stoi(clean, &consumed);
    if (consumed != clean.size()) {
        throw std::invalid_argument("Invalid integer value for " + name + ".");
    }
    return value;
}

std::uint64_t parse_uint64_scalar(const std::string& text, const std::string& name) {
    const std::string clean = unquote(text);
    std::size_t consumed = 0;
    const std::uint64_t value = static_cast<std::uint64_t>(std::stoull(clean, &consumed));
    if (consumed != clean.size()) {
        throw std::invalid_argument("Invalid uint64 value for " + name + ".");
    }
    return value;
}

math::Vec3 parse_vec3_scalar(const std::string& text, const std::string& name) {
    std::string clean = unquote(text);
    clean = trim(clean);
    if (clean.size() < 2 || clean.front() != '[' || clean.back() != ']') {
        throw std::invalid_argument(name + " must be a bracketed 3-vector.");
    }
    clean = clean.substr(1, clean.size() - 2);

    std::vector<double> values;
    std::stringstream ss(clean);
    std::string item;
    while (std::getline(ss, item, ',')) {
        values.push_back(parse_double_scalar(trim(item), name));
    }
    if (values.size() != 3) {
        throw std::invalid_argument(name + " must contain exactly 3 values.");
    }
    return math::Vec3{values[0], values[1], values[2]};
}

Method parse_method(const Node& root) {
    const std::string method_text = unquote(require_scalar(require_key(root, "method", "root"), "method"));
    if (method_text == "wop") {
        return Method::Wop;
    }
    if (method_text == "wos") {
        return Method::Wos;
    }
    throw std::invalid_argument("method must be one of: wop, wos.");
}

FunctionKind parse_function_kind(const std::string& text) {
    const std::string clean = unquote(text);
    if (clean == "constant") {
        return FunctionKind::Constant;
    }
    if (clean == "x") {
        return FunctionKind::X;
    }
    if (clean == "y") {
        return FunctionKind::Y;
    }
    if (clean == "z") {
        return FunctionKind::Z;
    }
    if (clean == "coulomb") {
        return FunctionKind::Coulomb;
    }
    throw std::invalid_argument("function kind must be one of: constant, x, y, z, coulomb.");
}

solver::RMaxMode parse_r_max_mode(const std::string& text) {
    const std::string clean = unquote(text);
    if (clean == "escape") {
        return solver::RMaxMode::Escape;
    }
    if (clean == "project") {
        return solver::RMaxMode::Project;
    }
    throw std::invalid_argument("wop.r_max_mode must be one of: escape, project.");
}

GeometryConfig parse_geometry(const Node& root) {
    const Node& geometry_node = require_map(require_key(root, "geometry", "root"), "geometry");

    GeometryConfig geometry;
    geometry.interior_point = parse_vec3_scalar(
        require_scalar(require_key(geometry_node, "interior_point", "geometry"), "geometry.interior_point"),
        "geometry.interior_point");

    const Node& planes_node = require_sequence(require_key(geometry_node, "planes", "geometry"), "geometry.planes");
    if (planes_node.sequence.empty()) {
        throw std::invalid_argument("geometry.planes must not be empty.");
    }

    geometry.planes.reserve(planes_node.sequence.size());
    for (const Node& plane_node : planes_node.sequence) {
        const Node& plane_map = require_map(plane_node, "geometry.planes item");
        geometry.planes.push_back(PlaneConfig{
            parse_vec3_scalar(require_scalar(require_key(plane_map, "p", "geometry.planes"), "geometry.planes[].p"), "geometry.planes[].p"),
            parse_vec3_scalar(require_scalar(require_key(plane_map, "nu", "geometry.planes"), "geometry.planes[].nu"), "geometry.planes[].nu"),
        });
    }

    return geometry;
}

FunctionConfig parse_function_section(const Node& root, const std::string& section_name) {
    const Node& fn_node = require_map(require_key(root, section_name, "root"), section_name);
    FunctionConfig fn;
    fn.kind = parse_function_kind(require_scalar(require_key(fn_node, "kind", section_name), section_name + ".kind"));

    if (fn.kind == FunctionKind::Constant) {
        fn.value = parse_double_scalar(
            require_scalar(require_key(fn_node, "value", section_name), section_name + ".value"),
            section_name + ".value");
    } else if (fn.kind == FunctionKind::Coulomb) {
        fn.source = parse_vec3_scalar(
            require_scalar(require_key(fn_node, "source", section_name), section_name + ".source"),
            section_name + ".source");
    }

    return fn;
}

std::optional<double> parse_optional_radius(const std::string& text) {
    const std::string clean = unquote(text);
    if (clean == "null") {
        return std::nullopt;
    }
    const double value = parse_double_scalar(clean, "wop.r_max");
    return (value > 0.0) ? std::optional<double>(value) : std::nullopt;
}

WopConfig parse_wop_section(const Node& root) {
    const Node& wop_node = require_map(require_key(root, "wop", "root"), "wop");
    WopConfig cfg;
    cfg.r_max = parse_optional_radius(require_scalar(require_key(wop_node, "r_max", "wop"), "wop.r_max"));
    cfg.r_max_mode = parse_r_max_mode(require_scalar(require_key(wop_node, "r_max_mode", "wop"), "wop.r_max_mode"));
    cfg.r_max_factor = parse_double_scalar(require_scalar(require_key(wop_node, "r_max_factor", "wop"), "wop.r_max_factor"), "wop.r_max_factor");
    if (cfg.r_max_factor <= 1.0) {
        throw std::invalid_argument("wop.r_max_factor must be greater than 1.0.");
    }
    return cfg;
}

WosConfig parse_wos_section(const Node& root) {
    const Node& wos_node = require_map(require_key(root, "wos", "root"), "wos");
    WosConfig cfg;
    cfg.delta = parse_double_scalar(require_scalar(require_key(wos_node, "delta", "wos"), "wos.delta"), "wos.delta");
    cfg.rho_scale = parse_double_scalar(require_scalar(require_key(wos_node, "rho_scale", "wos"), "wos.rho_scale"), "wos.rho_scale");
    cfg.rho1_scale = parse_double_scalar(require_scalar(require_key(wos_node, "rho1_scale", "wos"), "wos.rho1_scale"), "wos.rho1_scale");
    if (cfg.delta <= 0.0) {
        throw std::invalid_argument("wos.delta must be positive.");
    }
    if (cfg.rho_scale <= 0.0) {
        throw std::invalid_argument("wos.rho_scale must be positive.");
    }
    if (cfg.rho1_scale <= 1.0) {
        throw std::invalid_argument("wos.rho1_scale must be greater than 1.0.");
    }
    return cfg;
}

std::string load_text_file(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::invalid_argument("Could not open config file: " + path.string());
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

}  // namespace

RuntimeConfig load_config_file(const std::filesystem::path& path) {
    const Node root = parse_yaml_subset(load_text_file(path));
    require_map(root, "root");

    RuntimeConfig cfg;
    cfg.method = parse_method(root);
    cfg.x0 = parse_vec3_scalar(require_scalar(require_key(root, "x0", "root"), "x0"), "x0");
    cfg.n = parse_int_scalar(require_scalar(require_key(root, "n", "root"), "n"), "n");
    cfg.seed = parse_uint64_scalar(require_scalar(require_key(root, "seed", "root"), "seed"), "seed");
    cfg.max_steps = parse_int_scalar(require_scalar(require_key(root, "max_steps", "root"), "max_steps"), "max_steps");
    if (const Node* u_inf_node = find_key(root, "u_inf"); u_inf_node != nullptr) {
        cfg.u_inf = parse_double_scalar(require_scalar(*u_inf_node, "u_inf"), "u_inf");
    }

    if (cfg.n <= 0) {
        throw std::invalid_argument("n must be positive.");
    }
    if (cfg.max_steps <= 0) {
        throw std::invalid_argument("max_steps must be positive.");
    }

    cfg.geometry = parse_geometry(root);

    cfg.boundary = parse_function_section(root, "boundary");

    if (const Node* reference_node = find_key(root, "reference"); reference_node != nullptr) {
        require_map(*reference_node, "reference");
        cfg.reference = parse_function_section(root, "reference");
    }

    if (cfg.method == Method::Wop) {
        cfg.wop = parse_wop_section(root);
    } else {
        cfg.wos = parse_wos_section(root);
    }

    return cfg;
}

}  // namespace wop::config
