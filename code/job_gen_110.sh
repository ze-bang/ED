source ./script_venv/bin/activate
mkdir h110
cd h110
mkdir Jpm_0.0
cd Jpm_0.0
for i in {1..50}
do
    mkdir h_$i
    python3 ../../helper.py -0.001 1.0 -0.001 $i 0 3 50 1 1 0 h_$i
    cd h_$i
    ../../../../HPhi.build/src/HPhi -e ./namelist_cg.def
    cd ..
    python3 ../../read.py h_$i 110
done
cd ..
python3 ../read_energy.py Jpm_0.0 50 0 3 0.0

