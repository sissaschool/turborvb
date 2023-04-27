In order to translate a subroutine modify the three subroutin below in /src:
makefun_bumb.f90   
makefun.f90       
makefun_pbc.f90  

Execute in the present directory:
command_bump.sh   for         makefun_bumb.f90                    
command.sh        for         makefun.f90 
command_pbc.sh    for         makefun_pbc.f90                           

the three scripts will copy the corresponding 
   makefun*_b.f90 
files into the directory:
   trunk/AD/reverse_cell/

After that execute

./make0branch 

This will create makefun*0 and correspondin makefun*0_b for length 1 loop 

and TurboRVB  is ready to work for computing forces.

THE PART BELOW IS NECESSARY ONLY IF SOMETHING WENT WRONG.
You should check by hand if the actions below were correctly handled by the 
scripts.

#######################################
The following part is automatically executed by the scripts.
Be careful using command_pbc: 
the input variables 
   cellscale,sinphase,cosphase,rphase
are added, in the subroutine arguments, after:
   ,c,rmucos,rmusin)
If these variables will be changed, the script command_pbc
will not work correctly.
#######################################



NB for makefun_pbc.f90  

There is a little more work before execution in this case:

Remove the lines 5-6, namely:
         use Cell,      only:cellscale 
         use mod_twisted,       only:sinphase,cosphase,rphase
  
add   as input variables:

      cellscale,sinphase,cosphase,rphase) 

and define them with the two lines:

         double precision, intent(in) :: cellscale(3),rphase(3)  
         double precision :: sinphase(3,nion,0:indtm),cosphase(nion,0:indtm)

END NB for makefun_pbc.f90



After execution in the directory /output:

Edit the corresponding make*_b.f90 file and:

Replace : USE NAMEMODULE_B  --->  USE NAMEMODULE
e.g. (USE CONSTANTS_B     -->  USE CONSTANTS).

Then edit the *_b.f90 file and add the following two lines:

  integer :: adi4ibuf,adr8ibuf,adi4buf(1024)
  real*8  :: adr8buf(1024)


Also before the line:

SELECT CASE (iopt)

add two lines:

  adi4ibuf=1
  adr8ibuf=1


Finally  copy  the modified *_b.f90 file from the 
AD/compmake/output directory to the  AD/reverse_cell directory.


At the en

Then it should work.




