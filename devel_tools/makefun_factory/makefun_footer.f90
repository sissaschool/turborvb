case default
write(6,*) 'WARNING makefun: orbital',iopt,'not found'

iflagerr=1

end select
! ** ** ** ** ** ** ** END OF JASTROW ORBITALS ** ** ** ** ** ** ** ** *

return
end
