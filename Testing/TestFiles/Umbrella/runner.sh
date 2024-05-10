#!/bin/bash
echo "Submitting equil_1"
cp sub.sh run.sh
sed -i "s/NAME/equil_1/g" run.sh
ID=$(sbatch run.sh | awk {print $1}'')
echo "equil_1 ID is $ID"
echo "equil_1 ID is $ID" >> SLURMID.dat

    echo  "Submitting equil_2"
cp sub.sh run.sh
sed -i "s/NAME/equil_2/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}'))
echo "equil_2 ID is $ID"
echo "equil_2 ID is $ID" >> SLURMID.dat

echo  "Submitting equil_3"
cp sub.sh run.sh
sed -i "s/NAME/equil_3/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}'))
echo "equil_3 ID is $ID"
echo "equil_3 ID is $ID" >> SLURMID.dat

echo  "Submitting equil_4"
cp sub.sh run.sh
sed -i "s/NAME/equil_4/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}'))
echo "equil_4 ID is $ID"
echo "equil_4 ID is $ID" >> SLURMID.dat

echo  "Submitting equil_5"
cp sub.sh run.sh
sed -i "s/NAME/equil_5/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}'))
echo "equil_5 ID is $ID"
echo "equil_5 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_1"
cp sub.sh run.sh
sed -i "s/NAME/prod_1/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_1 ID is $ID"
echo "prod_1 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_2"
cp sub.sh run.sh
sed -i "s/NAME/prod_2/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_2 ID is $ID"
echo "prod_2 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_3"
cp sub.sh run.sh
sed -i "s/NAME/prod_3/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_3 ID is $ID"
echo "prod_3 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_4"
cp sub.sh run.sh
sed -i "s/NAME/prod_4/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_4 ID is $ID"
echo "prod_4 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_5"
cp sub.sh run.sh
sed -i "s/NAME/prod_5/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_5 ID is $ID"
echo "prod_5 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_6"
cp sub.sh run.sh
sed -i "s/NAME/prod_6/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_6 ID is $ID"
echo "prod_6 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_7"
cp sub.sh run.sh
sed -i "s/NAME/prod_7/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_7 ID is $ID"
echo "prod_7 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_8"
cp sub.sh run.sh
sed -i "s/NAME/prod_8/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_8 ID is $ID"
echo "prod_8 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_9"
cp sub.sh run.sh
sed -i "s/NAME/prod_9/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_9 ID is $ID"
echo "prod_9 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_10"
cp sub.sh run.sh
sed -i "s/NAME/prod_10/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_10 ID is $ID"
echo "prod_10 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_11"
cp sub.sh run.sh
sed -i "s/NAME/prod_11/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_11 ID is $ID"
echo "prod_11 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_12"
cp sub.sh run.sh
sed -i "s/NAME/prod_12/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_12 ID is $ID"
echo "prod_12 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_13"
cp sub.sh run.sh
sed -i "s/NAME/prod_13/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_13 ID is $ID"
echo "prod_13 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_14"
cp sub.sh run.sh
sed -i "s/NAME/prod_14/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_14 ID is $ID"
echo "prod_14 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_15"
cp sub.sh run.sh
sed -i "s/NAME/prod_15/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_15 ID is $ID"
echo "prod_15 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_16"
cp sub.sh run.sh
sed -i "s/NAME/prod_16/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_16 ID is $ID"
echo "prod_16 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_17"
cp sub.sh run.sh
sed -i "s/NAME/prod_17/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_17 ID is $ID"
echo "prod_17 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_18"
cp sub.sh run.sh
sed -i "s/NAME/prod_18/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_18 ID is $ID"
echo "prod_18 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_19"
cp sub.sh run.sh
sed -i "s/NAME/prod_19/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_19 ID is $ID"
echo "prod_19 ID is $ID" >> SLURMID.dat

echo  "Submitting prod_20"
cp sub.sh run.sh
sed -i "s/NAME/prod_20/g" run.sh
sed -i "s/#dep/#SBATCH --dependency=afterok:$ID/g" run.sh
ID=$(sbatch run.sh | awk '{print $1}')
echo "prod_20 ID is $ID"
echo "prod_20 ID is $ID" >> SLURMID.dat


