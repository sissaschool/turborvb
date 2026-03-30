subroutine print_version
    use logger_io, only: log_info
    implicit none
    character(LEN=6) :: version_number = '1.0.0'
    character(LEN=40) :: git_revision = 'unknown'

    call log_info("----------------------------------------------------------------------------")
    call log_info("TurboRVB version ", version_number, " git rev. ", git_revision)
    call log_info("")
    call log_info("  Ab-initio Quantum Monte Carlo Package")
    call log_info("")
    call log_info("    Developer: Sandro Sorella")
    call log_info("    Website: https://turborvb.sissa.it")
    call log_info("    GitHub: https://github.com/sissaschool/turborvb")
    call log_info("    Project PIs: Michele Casula and Kosuke Nakano")
    call log_info("    Contacts: michele.casula@gmail.com and kousuke_1123@icloud.com")
    call log_info("")
    call log_info("  When you publish a paper using TurboRVB, please cite the following paper.")
    call log_info("")
    call log_info("    TurboRVB: a many-body toolkit for ab initio electronic simulations,")
    call log_info("    K. Nakano*, C. Attaccalite, M. Barborini, L. Capriotti, M. Casula*,")
    call log_info("    E. Coccia, M. Dagrada, Y. Luo, G. Mazzola, A. Zen, and S. Sorella*,")
    call log_info("    J. Chem. Phys. 152, 204121 (2020), doi:10.1063/5.0005037")
    call log_info("")
    call log_info("----------------------------------------------------------------------------")
end subroutine print_version
