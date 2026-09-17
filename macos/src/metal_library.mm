/*
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "metal_library.h"
#include "metal_library_src.h"

#import <dispatch/dispatch.h>

id<MTLLibrary> OpenEMSMetalLibrary(id<MTLDevice> device, NSError** error)
{
	// The bytes are a static array, so no ownership transfer is needed.
	dispatch_data_t data = dispatch_data_create(metalLibrary, metalLibrary_len, NULL, ^{});
	return [device newLibraryWithData:data error:error];
}
