#pragma once

#include <memory>
#include <filesystem>

struct Config;
struct AMRState;

std::shared_ptr<AMRState> load_file(std::filesystem::path, Config const& c);

std::array<size_t, 3> mins(AMRState const&);
std::array<size_t, 3> maxs(AMRState const&);

bool value_at(AMRState const&, size_t i, size_t j, size_t k, float&);
