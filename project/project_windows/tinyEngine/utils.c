#include <stdio.h>
const char* getFilePath(const char* name, const char* type, char *path) {
	sprintf(path, "..\\..\\..\\resources\\image\\%s.%s", name, type);
	return path;
}