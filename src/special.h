#pragma once

struct Arguments;

/// Take a source volume, and generate meshes for all isovalues
/// the meshes are then merged. This is to handle the case where grid values are
/// actually identifiers, not some quantity
int all_iso_merge(Arguments& c);
