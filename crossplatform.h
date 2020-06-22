#ifndef CROSSPLATFORM_H
#define CROSSPLATFORM_H

#ifndef _WIN32

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#define mkdir(x) mkdir(x, 0777)

#else

#define _CRT_SECURE_NO_WARNINGS

#include <direct.h>

#define chdir _chdir
#define getcwd _getcwd
#define mkdir _mkdir



#endif // _WIN32

#endif // CROSSPLATFORM_H