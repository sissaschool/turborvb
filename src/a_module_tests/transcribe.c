/*
! Copyright (C) 2022 TurboRVB group
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include <stddef.h>
void transcribe_( double*A, int*b ) {
	void * ptr_before = (void*) A;
	void * ptr_after = NULL;

#pragma omp target data use_device_ptr(A)
	ptr_after = (void*) A;

	if (ptr_before != ptr_after) {
		*b = 0;
	}
	else {
		*b = 1;
	};
};
