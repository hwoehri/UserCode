#!/bin/sh

export STAGE_HOST=castorcms
export STAGE_SVCCLASS=cmst3

castorDir=/castor/cern.ch/user/h/hwoehri/cmst3/Polarization/Acceptance/JpsiGun_WithFSR/4Dec2010
cmsswDir=/afs/cern.ch/user/h/hwoehri/scratch0/Polarization/Acceptances/CMSSW_3_8_6/src

for ((job=1;job<10;job++));
  do
  echo "JOB "${job}
  name="jpsiGun_WithFSR_"${job}
#  displayfilename="display_"${name}".root"
  treefilename=${name}"_Tree.root"
  echo ${name}

  seed1=$(( ($job+1)*1307 ))
  sed -e "s/==SEED==/${seed1}/" jpsiGun_WithFSR_cfi_py_GEN.py > tmp_cfg

#Start to write the script
cat > job_${name}.sh << EOF
#!/bin/sh
cd $cmsswDir
eval \`scramv1 runtime -sh\`
cd -
#commande pour decoder le .cfg
cat > TEST_cfg.py << "EOF"
EOF

#Ajoute le .cfg au script
cat tmp_cfg>> job_${name}.sh

# On poursuit le script
echo "EOF" >> job_${name}.sh
cat >> job_${name}.sh << EOF
workDir=$PWD
cp ${cmsswdir}/jpsiGun_WithFSR_cfi_py_GEN.py .
echo "before running"
ls -l
cmsRun TEST_cfg.py >& ${name}.log
echo "after running"
ls -l

rfcp jpsiGun_Tree.root ${castorDir}/${treefilename}

rm -f jpsiGun_Tree.root

gzip ${name}.log
rfcp ${name}.log.gz $castorDir/log/
rm -f ${name}.log.gz

EOF
chmod 755 job_${name}.sh
bsub -q 8nh -J $name $PWD/job_${name}.sh

done

