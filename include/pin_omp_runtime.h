// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
//
// Pin the OpenMP runtime into the process on Windows.
//
// MATLAB unloads MEX files with FreeLibrary ("clear mex", recompilation, and
// process shutdown). When the last OpenMP-linked MEX unloads, the OpenMP
// runtime DLL's refcount hits zero and it is unmapped while its parked
// worker-pool threads are still alive; the next wake or teardown executes
// unmapped code and the process dies with 0xc0000005 (Windows Event Log:
// faulting module "VCOMP140.DLL_unloaded"). Pinning keeps the runtime
// resident for the process lifetime, so the worker pool stays valid across
// clear mex. GetModuleHandleEx is a no-op for runtimes that are not loaded.
#pragma once
#if defined(_WIN32) && !defined(NO_OMP)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
namespace gdsge_detail {
struct PinOmpRuntime {
    PinOmpRuntime() {
        HMODULE h = 0;
        GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_PIN, "vcomp140.dll", &h);   // MSVC
        GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_PIN, "libgomp-1.dll", &h);  // MinGW
    }
};
static PinOmpRuntime gdsge_pin_omp_runtime_instance;
}
#endif
