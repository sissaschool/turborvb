import glob

# Load header
with open('makefun_header.f90', 'r') as f:
    header = f.read()

# Load footer
with open('makefun_footer.f90', 'r') as f:
    footer = f.read()

# Load orbitals
# search for all file that fits name pattern orb_*.f90
orb_files = glob.glob('orb_*.f90')
# read files
orbitals = {f.replace("orb_", "").replace(".f90",""): open(f, 'r').read() for f in orb_files}

# Assemble output
with open('makefun_out.f90', 'w') as f:
    f.write(header)
    f.write('select case (iopt)\n')
    for k, v in orbitals.items():
        f.write(f'case ({k})\n')
        f.write(v)
    f.write(footer)

