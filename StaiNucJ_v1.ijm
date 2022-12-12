/* 
    StaiNucJ is an ImageJ macro developed to analyze different parameters
    in the cytoplasma and nuclei of cells.
    Copyright (C) 2022  Jorge Valero Gómez-Lobo.

    StaiNucJ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    StaiNucJ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//This macro has been developed by Dr Jorge Valero (jorgevalero@usal.es). 
//If you have any doubt about how to use it, please contact me.

//License
Dialog.create("GNU GPL License");
Dialog.addMessage("StaiNucJ  Copyright (C) 2022 Jorge Valero Gomez-Lobo.");
Dialog.setInsets(10, 20, 0);
Dialog.addMessage(" StaiNucJ comes with ABSOLUTELY NO WARRANTY; click on help button for details.");
Dialog.setInsets(0, 20, 0);
Dialog.addMessage("This is free software, and you are welcome to redistribute it under certain conditions; click on help button for details.");
Dialog.addHelp("http://www.gnu.org/licenses/gpl.html");
Dialog.show();








//*******************VARIABLES****************

//Variables of cells
sizecells=1000;
//Filter
gaussminCells=2;
gaussmaxCells=100;

//Nuclei variables
sizeNuc=500;
//Filter
closefilterRadius=4;
//Prominence for Find maxima
TolNuc=10;
//Size for refinement
sizeRefine=500;
//Band size for intensity measurement around the nuclei
bandsize=5;

//Measurements to be used
run("Set Measurements...", "area mean centroid integrated redirect=None decimal=2");

//Cleaning
run("Close All");
roiManager("reset");
waitForUser("Open the image");


//Image handling
name=getTitle();
selectWindow(name);
run("Make Composite");
Stack.setChannel(2);
run("Enhance Contrast...", "saturated=0.35");
Stack.setChannel(3);
run("Grays");
run("Enhance Contrast...", "saturated=0.35");

//Region ROI (first one, ROI 0)
setTool("polygon");
waitForUser("Draw the region of interest and press t");
roiManager("Add");
run("Select None");


//******************Cells processing********************
run("Duplicate...", "title=Cells duplicate channels=2");
run("Grays");
selectWindow("Cells");
//Difference of Gaussians filter
difGauss(gaussminCells, gaussmaxCells, "Cells");
//User thresholding assistance
run("Threshold...");
roiManager("Select", 0);
setAutoThreshold("Huang dark");
waitForUser("Set Threshold");
getThreshold(lowerCells, upperCells);

//**************Nuclei processing***********************
selectWindow(name);
run("Duplicate...", "title=Nuclei duplicate channels=3");
run("Grays");
run("Duplicate...", "title=Nucleiclose duplicate channels=3");
//Close filter for Nuclei
run("Minimum...", "radius="+closefilterRadius);
run("Maximum...", "radius="+closefilterRadius);
//User thresholding assistance
run("Threshold...");
roiManager("Select", 0);
setAutoThreshold("Mean dark");
waitForUser("Set Threshold for all nuclei");
getThreshold(lowerNuc, upperNuc);


//*******************Automatic identification of staining and nuclei**********
setBatchMode(true);

//**********************Staining identificaion****************************
//Staining ROIs Generation
selectWindow("DifGaussCells");
roiManager("Select", 0);
setThreshold(lowerCells, upperCells);
run("Analyze Particles...", "size="+sizecells+"-Infinity add");
RoisCells=roiManager("count");

//Closing images
selectWindow("DifGaussCells")
close();
selectWindow("Cells");
close();

//**********************Nuclei identificaion****************************
//Nuclei ROIs Generation (exclusively inside the staining)
selectWindow("Nucleiclose");
roiManager("Select", 0);
setThreshold(lowerNuc, upperNuc);
//Find maxima
run("Find Maxima...", "prominence="+TolNuc+" above output=[Segmented Particles]");
selectWindow("Nucleiclose Segmented");
roiManager("Select", Array.getSequence(RoisCells));
roiManager("XOR");
run("Make Inverse");
roiManager("Add");
roiManager("Select", newArray(0, RoisCells));
roiManager("AND");
setThreshold(1, 255);
//ROIS generation with Analyze particles
run("Analyze Particles...", "size="+sizeNuc+"-Infinity add");
roiManager("Deselect");
roiManager("Select", RoisCells);
roiManager("Delete");
RoisNuc=roiManager("count");
selectWindow("Nucleiclose");
close();

//Nuclei refinement (to better adjust to nuclei contour in original image)
selectWindow("Nuclei");
for (i=1; i<RoisCells; i++){
	roiManager("Select", i);
	setAutoThreshold("MinError dark");
	run("Analyze Particles...", "size="+sizeRefine+"-Infinity add");
}



//Deletion of not refined nuclei
selector=ArrayCreator(RoisCells, RoisNuc-1, 1);
roiManager("Select", selector);
roiManager("Delete");
RoisNuc=roiManager("count");

//Naming of Nuclei based on label
ROIrenamer(RoisCells, RoisNuc-1, "Nuc");
setBatchMode(false);
setBatchMode("show");




//*************CHANGE COLOR OF SELECTED NOT CONSIDERED ROIS*************
//Selection of Cell Markers

selector=ArrayCreator(RoisCells, RoisNuc-1, 1);
selectWindow("Nuclei");
roiManager("Deselect");
roiManager("Select", selector);
//Get info of nuclei
roiManager("Measure");
selectWindow("Results");
areaArr=Table.getColumn("Area");
meanArr=Table.getColumn("Mean");
Array.getStatistics(areaArr, minArea, maxArea, meanArea, stdDevArea);
Array.getStatistics(meanArr, minMean, maxMean, meanMean, stdDevMean);
selectWindow("Results");
run("Close");
do{
	selectWindow("Nuclei");
	roiManager("Select", Array.getSequence(roiManager("count")));
	roiManager("Set Color", "yellow"); 
	Plot.create("Results", "MeanIntensity", "Area");
	Plot.setLimits(0, maxMean+10, 0, maxArea+100);
	Plot.setColor("blue", "blue");
	Plot.add("circle", meanArr, areaArr);
	Plot.update;

	PlotSelector(false);
	selectWindow("Results");
	selectWindow("Nuclei");
	resetThreshold();
	roiManager("show all without labels");
	waitForUser("Check ROIs");
	repeat=getBoolean("Do you want to change limits?");
	
}  while(repeat==true);
PlotSelector(true);
RoisNuc=roiManager("count");
selectWindow("Results");
rename("Plot nuclei");

setBatchMode(true);
//*************** RESULTS TABLES *************************

//*******************Info nuclei and bands ******************

selectWindow(name);
run("Duplicate...", "title=Staining duplicate channels=2");

bandIntensities=newArray(RoisNuc-RoisCells);
NucArea=newArray(RoisNuc-RoisCells);
NucMean=newArray(RoisNuc-RoisCells);

selectWindow("Staining");
counter=0;
for (i=RoisCells; i<RoisNuc; i++){
	roiManager("Select", i); 
	getStatistics(areaNuc, meanNuc, min, max, std, histogram);
	NucArea[counter]=areaNuc;
	NucMean[counter]=meanNuc;
	run("Make Band...", "band="+bandsize);
	getStatistics(areaband, meanband, minband, maxband, std, histogram);
	bandIntensities[counter]=meanband;
	counter++;
}

Table.create("Nuclei and bands");
Table.setColumn("Area", NucArea);
Table.setColumn("Intensity Nuclei", NucMean);
Table.setColumn("Intensity Band", bandIntensities);

selector=ArrayCreator(RoisCells, RoisNuc-1, 1);
roiManager("Select", selector);
roiManager("Measure");
selectWindow("Results");
X=Table.getColumn("X");
Y=Table.getColumn("Y");
selectWindow("Results");
run("Close");


selector=ArrayCreator(1, RoisCells-1, 1);
roiManager("Select", selector);
roiManager("Measure");

AreaSt=Table.getColumn("Area");
IntensitySt=Table.getColumn("Mean");
IntegratedDn=Table.getColumn("IntDen");
RawIntegratedDn=Table.getColumn("RawIntDen");

numbernuc=newArray(RoisCells-1);

for (i=1; i<RoisCells; i++){
	counter=0;
	nucleicount=0;
	for (ii=RoisCells; ii<RoisNuc; ii++){
		roiManager("Select", i);
		if (selectionContains(parseFloat(X[counter]), parseFloat(Y[counter]))) nucleicount++;
		counter++;
	}
	numbernuc[i-1]=nucleicount;
}

Table.create("Cells/clusters staining and nuclei");
Table.setColumn("Area", AreaSt);
Table.setColumn("Intensity", IntensitySt);
Table.setColumn("Integrated Density", IntegratedDn);
Table.setColumn("Raw Integrated Density", RawIntegratedDn);
Table.setColumn("Nuclei", numbernuc);

selectWindow("Results");
run("Close");

roiManager("Select", 0);
roiManager("Measure");
selectWindow("Results");
rename("Results Selected Region");

setBatchMode(false);
setBatchMode("show");

function difGauss(min, max, extra){
	
	run("Duplicate...", "title=min duplicate");
	run("Gaussian Blur...", "sigma="+min+" stack");
	run("Duplicate...", "title=max duplicate");
	run("Gaussian Blur...", "sigma="+max+" stack");

	imageCalculator("Subtract stack", "min","max");
	selectWindow("max");
	close();
	selectWindow("min");
	rename("DifGauss"+extra);	
}


function ArrayCreator(from, to, step){
	if ((to>=from && step>=0) || (to<from && step<0)){
		if (to>=from && step>=0) {
			Arr=newArray();
			counter=0;
			for (i=from; i<=to; i=i+step){
				Arr[counter]=i;
				counter++;
			}
			return Arr;
		}
		if (to<from && step<0){
			Arr=newArray();
			counter=0;
			for (i=from; i>=to; i=i+step){
				Arr[counter]=i;
				counter++;
			}
			return Arr;
		} 
	}
	else return newArray("To > from and step <0 or To<from and step>0");	
}


function ROIrenamer(ini, fin, prefix){
	for (i=ini; i<=fin; i++){
		roiManager("Select", i);
		number=roiManager("index")+1;
		roiManager("Rename", prefix+number);
	}
	
}

function PlotSelector(erase){
	//ROIs exclusion using plot
	selectWindow("Results");
	if (erase==false){
		do {
			waitForUser("Create ROI for qualifying some dots, then click OK (press shift to add more than one ROI)");
		} while(selectionType == -1)
	}
	
	qualifiedX = newArray(meanArr.length);
	qualifiedY = newArray(areaArr.length);
	roiTochange = newArray(areaArr.length);
	q = 0;
	   for (jj = 0; jj <meanArr.length; jj++){
	      x = meanArr[jj];
	      y = areaArr[jj];
	      toUnscaled(x, y);
	      if (selectionContains(x, y)){
	         qualifiedX[q] = meanArr[jj];
	         qualifiedY[q] = areaArr[jj];
	         roiTochange[q] = jj+RoisCells;
	         q++;
	      }
	   }
	   qualifiedX = Array.trim(qualifiedX, q);
	   qualifiedY = Array.trim(qualifiedY, q);
	   roiTochange = Array.trim(roiTochange, q);
	
	   Plot.create("Results", "Mean Intensity", "Area");
	   Plot.setLimits(0, maxMean+10, 0, maxArea+100);
	   Plot.setColor("blue", "blue");
	   Plot.add("circles", meanArr, areaArr); 
	   Plot.setColor("red");
	   Plot.add("circles", qualifiedX, qualifiedY); 
	   Plot.update;
	   if (erase==false){
	   	selectWindow("Nuclei");
	   	roiManager("Select", roiTochange);
	   	roiManager("Set Color", "red");  
	   }
	   if (erase==true){
	   	selectWindow("Nuclei");
	   	roiManager("Deselect");
	   	roiManager("Select", roiTochange);
	   	roiManager("Delete");
	   	roiManager("show all without labels");
	   }
	   
}
