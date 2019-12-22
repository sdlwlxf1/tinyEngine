#include <stdio.h>
const char* getFilePath(const char* name, const char* type) {
	static char path[100];
	if(strcmp(type, "png") == 0)
		sprintf_s(path, sizeof(path), "../../../resources/image/%s.%s", name, type);
	else if(strcmp(type, "obj") == 0)
		sprintf_s(path, sizeof(path), "../../../resources/obj/%s.%s", name, type);
	return path;
}