subroutine print_version
    implicit none
    character(LEN=6) :: version_number = '0.0.1'
    character(LEN=40) :: git_revision = 'unknown'

    write (6, *) "----------------------------------------------------------------------------"
    write (6, *) "TurboRVB version ", version_number, " git rev. ", git_revision
    write (6, *) ""
    write (6, *) "  Ab-initio Quantum Monte Carlo Package"
    write (6, *) ""
    write (6, *) "    Developer: Sandro Sorella"
    write (6, *) "    Website: https://turborvb.sissa.it"
    write (6, *) "    GitHub: https://github.com/sissaschool/turborvb"
    write (6, *) "    Project PIs: Michele Casula and Kosuke Nakano"
    write (6, *) "    Contacts: michele.casula@gmail.com and kousuke_1123@icloud.com"
    write (6, *) ""
    write (6, *) "  When you publish a paper using TurboRVB, please cite the following paper."
    write (6, *) ""
    write (6, *) "    TurboRVB: a many-body toolkit for ab initio electronic simulations,"
    write (6, *) "    K. Nakano*, C. Attaccalite, M. Barborini, L. Capriotti, M. Casula*,"
    write (6, *) "    E. Coccia, M. Dagrada, Y. Luo, G. Mazzola, A. Zen, and S. Sorella*,"
    write (6, *) "    J. Chem. Phys. 152, 204121 (2020), doi:10.1063/5.0005037"
    write (6, *) ""
    write (6, *) "----------------------------------------------------------------------------"
end subroutine print_version
