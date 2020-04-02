
#if defined(__unix__) || defined(__linux__)

#include <pthread.h>

void pin_to_core(size_t core)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
}

#endif


#if defined(_WIN32) || defined(_WIN64)
#define NOMINMAX
#include <Windows.h>

//TODO: Test this implementation
void pin_to_core(size_t core)
{
	DWORD_PTR mask = 1;
	mask <<= core;
	SetThreadAffinityMask(GetCurrentThread(), mask);
}
#endif
