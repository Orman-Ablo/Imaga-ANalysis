//Setting up channels and duplicates to preserve the original images
//This makes a max-intensity image from the screening image
roiManager("reset");
rename("Cells");
run("Z Project...", "projection=[Max Intensity]");
//This splits up the fluorescent markers into Cadherin, DAPI and VWF
run("Split Channels");
selectImage("C1-MAX_Cells");
run("Duplicate...", "title=VWF");
selectImage("C2-MAX_Cells");
run("Duplicate...", "title=Cadherin");
selectImage("C3-MAX_Cells");
run("Duplicate...", "title=DAPI");
//Creating masks from the nuclei: to perform marker-based watershed
setOption("BlackBackground", true);
run("Convert to Mask");
run("Analyze Particles...", "size=500-Infinity pixel display add");
roiManager("Combine");
//these next few lines removes background noise, and fills out some selections
//This will help the watershed work smoother
setBackgroundColor(0, 0, 0);
run("Clear Outside");
setForegroundColor(255, 255, 255);
run("Fill", "slice");
selectImage("Cadherin");
// A gaussian blur helps remove background noise
//And identify boundaries.
run("Gaussian Blur...", "sigma=2");
//Marker-controlled watershed takes your image and a marker image 
//it uses the marker image as a collection of 'seed points' from where it can 
// create a watershed onto the chosen image. This allows for effective membrane 
//membrane identification
run("Marker-controlled Watershed", "input=Cadherin marker=DAPI mask=None compactness=0 binary calculate use");
//We threshold the watershed to only exclude the boundaries(intensity 0). 
//This provides a neutral 'map' of the cellular membranes
run("Duplicate...", "title=[Cell Membranes]");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
run("Convert to Mask");
//Select all nuclear selections
count = roiManager("count");
//count upp all selections (Nuclei) in the roiManager
//(We could also save them if desired)
array = newArray(count);
  for (i=0; i<array.length; i++) {
      array[i] = i;
  }
//Delete them: we'll want ot cycle over the membrane selections later
roiManager("select", array);
roiManager("Delete");
run("Close");
//Make membrane selections of only the full cells: edges excluded
//This prevents us from including only part of a cell
//(If we include, Shouldn't the reduced area correspond to a lower amount of bodies too: 
//It should work out?)
selectImage("Cell Membranes");
setThreshold(1, 255);
setOption("BlackBackground", true);
run("Convert to Mask");
//add or remove 'exclude' to choose whether you want to look at all cells
//or only whole cells
run("Analyze Particles...", "size=500-Infinity pixel display exclude add");
//Just to collect Cell Info, we need to clear the results Table
Table.create("Results");
//Measure the prescence of VWF in the full image
selectImage("VWF");
roiManager("Select", 0);
roiManager("Add");
roiManager("deselect");
roiManager("Measure");
nselections = roiManager("count");
setAutoThreshold("MaxEntropy dark");
run("Set Measurements...", "area mean standard perimeter shape feret's integrated median area_fraction display redirect=None decimal=3");
//This for loop analyzes the particles within each selection (cell)
for (i = 0; i < nselections; i++) {
	selectImage("VWF");
	roiManager("select", i);
	roiManager("Rename", "Cell "+(i+1));
	
	run("Analyze Particles...", "size=3-200 pixel display exclude add");

//Unused code below: previously used to make duplicate images per cell
//run("Duplicate...", "title=Cell");
//rename("Cell "+(i+1));
//run("Clear Outside");
}
//Closing some unneccesary images
selectImage("C1-MAX_Cells");
close();
selectImage("C2-MAX_Cells");
close();
selectImage("C3-MAX_Cells");
close();
selectImage("DAPI");
close();
selectImage("Cadherin");
close();
selectImage("Cell Membranes");
close();

//ignore below:

//For each cell, we will measure the longest length of objects using 'ferret's diameter'
//roiManager("deselect");
//for (i = 0; i < nselections; i++) {
//Table.create("VWF info cell " +(i+1));
//selectImage("Cell "+(i+1));
//roiManager("select", i);
//setAutoThreshold("MaxEntropy dark");
//run("Analyze Particles...", "size=1-200 pixel display exclude add");
//close();
//}

