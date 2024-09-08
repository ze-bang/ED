source ./script_venv/bin/activate
mkdir h111
cd h111
mkdir Jpm_0.025
cd Jpm_0.025
for i in {1..50}
do
    mkdir h_$i
    python3 ../../helper.py -0.05 1.0 -0.05 $i 0 1 50 1 1 1 h_$i
    cd h_$i
    ../../../../HPhi.build/src/HPhi -e ./namelist_cg.def
    cd ..
    python3 ../../read.py h_$i 111
done
cd ..
python3 ../read_energy.py Jpm_0.025 50 0 1 0.025

