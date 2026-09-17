/*
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef METAL_LIBRARY_H
#define METAL_LIBRARY_H

#import <Metal/Metal.h>

//! The precompiled kernel library embedded in libopenEMS, or nil with \a error
//! set. Shaders are compiled at build time; the runtime never compiles.
id<MTLLibrary> OpenEMSMetalLibrary(id<MTLDevice> device, NSError** error);

#endif // METAL_LIBRARY_H
