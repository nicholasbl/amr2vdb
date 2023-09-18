#pragma once

struct Arguments;

/// Turn an AMR grid into a VDB volume
int amr_to_volume(Arguments const& c);
