#! /bin/bash
# DTI preprocessing with AAL 116 or Desikan-Killiany 68 parcellation
# batch processing
## ("#!" is an operator called shebang which directs the script to the interpreter location). 
## Assuming we are in the required directory which has all the dti files and T1.w file with this bash script file##
##To run this file, open terminal and type "chmod +x run_dti.sh" to make this file executable and run this file , type ./run_dti.sh##
# To run a single subject, arrange following data in a folder: dti data = data.nii.gz, bvals, bvecs, T1 and ROI-atlas (AAL)

for j in ANAL ;do # change 'sub-*' based on on your subject wise folder name
echo $j
cd $j
d=data.nii.gz*
bvc=bvecs
bvl=bval*
t=T1w.nii.gz*
#dap=AP.nii.gz*
#dpa=PA.nii.gz*
#### Step 1 - Preprocessing of DTI data ######

### Step 1.1 - Conversion of DTI file formats into mif (MRtrix Image Format)
mrconvert $d sub.mif -fslgrad $bvc $bvl -force
## few  seconds

#mrview sub${i}.mif
### Step 1.2 - Denoising the data
dwidenoise sub.mif sub_den.mif -noise noise.mif -force
## 2 mins

##Step 1.3 -distortion map formation

#mrconvert $dap b0_AP.mif
#mrconvert $dpa b0_PA.mif
#mrcat b0_AP.mif b0_PA.mif -axis 3 b0_pair.mif -force 

### Step 1.4 - Motion and distortion correction
dwifslpreproc sub_den.mif sub_preproc.mif -rpe_none -pe_dir ap -readout_time 0.038

### Step 1.5 - Bias Field correction
dwibiascorrect ants sub_preproc.mif sub_unbiased.mif -force
## 1 min

### Step 1.6 - Brain mask estimation 
#mrconvert sub_unbiased.mif sub_unbiased.nii -export_grad_fsl bvecs bvals -force
#bet2 sub_unbiased.nii sub_masked.nii -m -f 0.3 # threshold needs to be adjusted based on your image
mri_synthstrip  -i T1.nii.gz -o T1_brain.nii.gz
5ttgen fsl T1.nii.gz 5tt.mif  #Five Tissue Types (5TT) segmentation

#### Step 2- T1w registration in DTI space
#Brain mask estimation on b0 space (for that first create the b0 file from DTI first volume)
fslroi data.nii.gz b0.nii.gz 0 1
bet b0.nii.gz b0_brain -f 0.25 -R -g 0.1
fslmaths b0_brain.nii.gz -bin b0_brain_mask
# T1w registration in DTI space using fsl-flirt
#flirt -in T1_brain.nii.gz -ref sub_masked.nii -omat t12b0.mat -out t12b0.nii.gz -dof 6
#flirt -in T1w_brainmask.nii.gz -ref sub_masked.nii.gz -applyxfm -init t12b0.mat -interp nearestneighbour -out mask_on_b0.nii.gz
flirt -ref T1_brain.nii.gz -in b0_brain.nii.gz -omat b02t1.mat -dof 6 -interp nearestneighbour
transformconvert b02t1.mat b0_brain.nii.gz T1_brain.nii.gz flirt_import diff2struct_mrtrix.txt
convert_xfm -omat t12b0.mat -inverse b02t1.mat
mrtransform 5tt.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif
5tt2gmwmi 5tt_coreg.mif gmwmseed_coreg.mif
run_first_all -i T1_brain.nii.gz -o T1_first -b 

#### Step 3- Fiber Orientation Distribution (FOD)
### Step 3.1 - Response function estimation
dwi2response dhollander  sub_unbiased.mif wm.txt gm.txt csf.txt -mask b0_brain_mask.nii.gz -voxels voxels.mif -force
## 30-45min
### Step 3.2 - Estimation of Fiber Orientation Distribution (FOD)
dwi2fod msmt_csd sub_unbiased.mif  -mask b0_brain_mask.nii.gz wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif -force
## 6 min
### Step 3.3 - Intensity Normalization 
#mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif -mask b0_brain_mask.nii.gz -force

### Step 4- Creating streamlines
tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmseed_coreg.mif  -maxlength 250 -cutoff 0.06 -select 10000000 wmfod.mif tracks_10M.tck -force
## 45min
### Step 4.1 - Reducing the number of streamlines

tcksift2 -act 5tt_coreg.mif -out_coeffs sift_coeffs.txt -nthreads 8 tracks_10M.tck wmfod.mif sift_1M.txt -force # [17 minutes]

#### Step 5- MNI template (which also contains non-brain tissue) to T1w space
#flirt -in MNI152_T1_1mm.nii.gz -ref T1w_preproc.nii.gz -omat mni2t1.mat -out mni2t1_affine.nii.gz -dof 12
flirt -in aal.nii -ref T1_brain.nii.gz -omat aal2t1.mat -out aal2t1_affine.nii.gz -dof 12 -interp nearestneighbour

#fnirt --in=mni2t1_affine.nii.gz --ref=T1w_preproc.nii.gz --cout=mni2t1_transf --refmask=T1w_brainmask.nii.gz  --iout=mni2t1_fnirted.nii.gz

##### Step 6- Transform the artlas (e.g., aal, Shen, schaefer) parcellation to the T1w space
#applywarp --ref=T1w_preproc.nii.gz --in= Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz --warp=mni2t1_transf.nii.gz  --premat=mni2t1.mat --out=Schaefer400_on_T1.nii.gz --interp=nn

##### Step 7-  Transform the Schafer parcel to the DWI_b0 space
flirt -in aal2t1_affine.nii.gz -ref b0_brain.nii.gz -applyxfm -init t12b0.mat -interp nearestneighbour -out aal2t1_b0.nii.gz

### Step 8- Generating connectome (ROI-to-ROI connectivity matrix)
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_1M.txt tracks_10M.tck aal2t1_b0.nii.gz SC.csv -out_assignment assignments.csv -force

### Step 9- Sreate specific region (e.g., Thalamus) tractography
fslsplit T1_first_all_fast_origsegs.nii.gz
fslmaths vol0006.nii.gz -add vol0013.nii.gz -bin thalamus_mask # specifay spesified ROIs .nii file here
#fsleyes thalamus_mask.nii.gz 
flirt -in thalamus_mask.nii.gz -ref b0_brain.nii.gz -applyxfm -init t12b0.mat -interp nearestneighbour -out thal2b0.nii.gz 
tckedit -include thal2b0.nii.gz tracks_10M.tck thalamus.tck

### View the results in MRview or TrackVis or DTIstudio
cd ..
done






