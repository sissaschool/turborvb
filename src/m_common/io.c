/*
Copyright (C) 2022 TurboRVB group based on code by
Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

*/

#if defined RISC
#define F90_FUNC_(name, NAME) name
#else
#define F90_FUNC_(name, NAME) name##_
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

void F90_FUNC_(imkdir, IMKDIR)(char *name)
{
        struct stat buf;
        if (!*name || stat(name, &buf) == 0)
                return;
        mkdir(name, 0775);
}

void F90_FUNC_(ichdir, ICHDIR)(char *name)
{
        chdir(name);
}

void F90_FUNC_(irename, IRENAME)(char *namein, char *nameout)
{
        rename(namein, nameout);
}

void F90_FUNC_(iremove, IREMOVE)(char *name)
{
        remove(name);
}

void F90_FUNC_(isystem, ISSYSTEM)(char *name, int *ierr)
{
        *ierr = system(name);
}

void F90_FUNC_(igetcwd, IGETCWD)(char *name, int *ln)
{
        getcwd(name, 256);
        *ln = strlen(name);
}

void F90_FUNC_(igethname, IGETHNAME)(char *name, int *ln)
{
        gethostname(name, 256);
        *ln = strlen(name);
}
