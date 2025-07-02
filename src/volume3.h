#pragma once

struct Arguments;

/// Turn an AMR grid into a series of VDB Volumes, one for each level
int amr_to_volume_sets(Arguments const& c);
