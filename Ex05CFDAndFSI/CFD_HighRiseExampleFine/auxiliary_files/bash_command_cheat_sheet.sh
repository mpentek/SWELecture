
# access to the computers
# stud1
ssh wise2021@129.187.141.198 
# or  
ssh wise2021@stud1.st.bv.tum.de

# stud11
ssh wise2021@129.187.141.111 
# or  
ssh wise2021@stud11.st.bv.tum.de

# scp recursive copy
# of the folder cfd.gid of your local machine 
# to the folder /home/WiSe2021/Documents/
# on the remote machine wise2021@stud1.st.bv.tum.de
scp -r cfd.gid/ wise2021@stud1.st.bv.tum.de:/home/WiSe2021/Documents/
# to check content
# login via ssh to wise2021@stud1.st.bv.tum.de
ssh wise2021@stud1.st.bv.tum.de
# navigate to folder on the remote machine and check content
cd Documents/cfd.gid/
# list content
ls
# structured with more detail
ll -l
# on the remote machine test in a screen and check the simul.log in VS Code
# screen executes the series of bash commands in the string 
# already in a detached state
screen -dm bash -c 'python3 MainKratosCustom.py | tee simul.log'
killall screen python3
# or in the terminal but here there is the danger 
# that some zombie processes might remain in the memory after killing it Ctrl + Z
python3 MainKratosCustom.py
killall python3

# check forces
cd results/ascii_output/forces/
python3 check_forces.py
python3 convert_kratos_to_paroptbeam.py 

# add cp data to h5 of deck results
cd ../../
python3 add_cp_to_h5.py

# create xdmf for h5 to be viewed with Paraview
cd hdf5_output/domain/
convertH5toXdmf fluid_computational_model_part-0.00.h5

cd ../deck/
convertH5toXdmf NoSlip3D_structure_T-0.00.h5

# start the final run and check with htop or top if processor usage increases
screen -dm bash -c '{ time taskset -c 1-6 python3 MainKratosCustom.py | tee simul.log ; } 2> time.log ; { tail -n 25 simul.log && cat time.log ; } | mail -s "Simulation finished on $(hostname)" someone@tum.de '
#top
htop

# when done with testing or simulation copy data to your local machine
# in the current working directory
scp -r WiSe2021@stud1.st.bv.tum.de:/home/WiSe2021/Documents/cfd.gid/ .