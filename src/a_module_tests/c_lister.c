#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include "c_lister.h"

extern void list_directory(char *path, int path_length, char *files, int *num_files, char *prefix, int prefix_length) {

    struct dirent *entry;

    char cpath[path_length + 1];
    strncpy(cpath, path, path_length);
    cpath[path_length] = '\0';

    DIR *dir = opendir(cpath);

    int count = 0;

    // Clean files
    int i;
    for (i = 0; i < MAX_FILES * MAX_FILE_LENGTH; i++) {
	files[i] = ' ';
    }

    while ((entry = readdir(dir)) != NULL && count < MAX_FILES) {
        if (entry->d_name[0] == '.') {
            continue;
        }
        if (strlen(entry->d_name) < prefix_length) {
            continue;
        }
        bool match = true;
        for (i = 0; i < prefix_length; i++) {
            if (entry->d_name[i] != prefix[i]) {
        	match = false;
        	break;
            }
        }
        if (!match) {
            continue;
        }
        //if (entry->d_name+strlen(entry->d_name) - 1 != "t") continue;
        //if (entry->d_name+strlen(entry->d_name) - 2 != "a") continue;
        //if (entry->d_name+strlen(entry->d_name) - 3 != "d") continue;
        //if (entry->d_name+strlen(entry->d_name) - 4 != ".") continue;
        
        strcpy(files + count * MAX_FILE_LENGTH, entry->d_name);
        count++;
    }
    *num_files = count;
    
    closedir(dir);
}




   
