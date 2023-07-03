#!/usr/bin/env python3

# Copyright (C) 2022 TurboRVB group
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Kosuke Nakano created on 21st Feb. 2019.

# program:read_colums.py
# purpose:read a column fron fort.21

# python modules
import pandas as pd
import numpy as np
import sys
import time


def main():
    # number_of_extracted_column=int(sys.argv[1])
    input_filename = str(sys.argv[1])
    output_filename = str(sys.argv[2])

    # reading fort.12
    print(
        "    [read_columns_fort21.py] Reading {} using pandas read_csv...".format(
            input_filename
        ),
        end="",
    )
    start = time.time()
    data_all = pd.read_csv(input_filename, header=None, delim_whitespace=True)
    print("done.")
    elapsed_time = time.time() - start
    print(
        "    [read_columns_fort21.py] The elapsed_time is {0}".format(
            elapsed_time
        )
        + "[sec]"
    )

    # rename columns
    columns = [i + 1 for i in range(len(data_all.columns) - 2)] + [
        "weight",
        "bin",
    ]
    data_all.columns = columns
    data_all["bin"] = data_all["bin"].astype(np.int64)
    # pd.options.display.float_format = '{:.12e}'.format
    # data_all.loc[:, [number_of_extracted_column, "weight", "bin"]].to_csv(output_filename, float_format='%.12e', header=False, index=None, sep=" ")

    # save all the columns
    print(
        "    [read_columns_fort21.py] Writing {} using pandas to_csv...".format(
            input_filename
        ),
        end="",
    )
    start = time.time()
    for num in [i + 1 for i in range(len(data_all.columns) - 2)]:
        data_all.loc[:, [num, "weight", "bin"]].to_csv(
            output_filename + "_{}".format(num),
            float_format="%.12e",
            header=False,
            index=None,
            sep=" ",
        )
    print("done.")
    elapsed_time = time.time() - start
    print(
        "    [read_columns_fort21.py] The elapsed_time is {0}".format(
            elapsed_time
        )
        + "[sec]"
    )


if __name__ == "__main__":
    main()
