//YggData is the ThompsonLab ImageJ macro for single cell analysis of IFA images

yeastSize="10-infinity";
yeastContrast="0.3";
yeastThreshold="102";

dapiExtract=true;
nucleusSize="5-infinity";
dnaMinThreshold="35";

conaExtract=true;
conaSize="10-infinity";
conaMinThreshold="40";
conaCirc = "0.1";


//	This macro first asks you to pick a nucleus picture to generate ROIs from
//	It then uses these ROIs to determine the raw NUCLEAR data from each image within the directory, and stores these CSVs in a newly created 'Nuclear' folder
//	Ygg will then back out and generate ROIs for each image independent of nuclear localization, and store these CSVs in a newly created 'WholeCell' folder
//	Following this macro, the imaGen() script to combine all the data into a new dataset in which each row is a single cell

// First, lets make sure the measurements are set for appropriate analyses. If changes are made here, they will not be included in the default imaGen() compilation

run("Set Measurements...", "area mean standard modal min centroid center perimeter feret's integrated median area_fraction display redirect=None decimal=9");

roi_image = File.openDialog("Choose a File");
input=getDirectory("current");
parent=File.getParent(input);
input=File.getParent(parent);

pList=getFileList(input);
for (i=0; i < pList.length; i++){
	if(!endsWith(pList[i], ".csv")){
		dList=getFileList(input+"/"+pList[i]);
		for (j=0; j < dList.length; j++){
		if(!endsWith(dList[j], ".csv")){
			pathway = input+"/"+pList[i]+"/"+dList[j];
			pathway2 = pathway+"/PNGS";
			iList = getFileList(pathway);
			if(!File.isDirectory(pathway+"/csvs")){
				File.makeDirectory(pathway+"/csvs");
			};
			inputc = pathway+"/csvs/";
//			open(pathway+"/bf.png");
//			run("Enhance Contrast...", "saturated="+yeastContrast);
//			run("8-bit");
		// Change this for deviations from default
//			setAutoThreshold("Default dark");
		//run("Threshold...");
//			setThreshold(0, yeastThreshold);
//			setOption("BlackBackground", false);
//			run("Convert to Mask");
//			run("Analyze Particles...", "size="+yeastSize+" add include exclude");
//			close();
//			open(pathway+"/bf.png");
//			run("8-bit");
//			roiManager("measure");
//			close();
//			saveAs("Results", inputc+"bf.csv");
//			run("Clear Results");
			
			if(conaExtract){
				roi_image = pathway2+"/cona.png";
				open(roi_image);
				run("8-bit");
				setAutoThreshold("Default dark");
				//run("Threshold...");
				setThreshold(conaMinThreshold, 255);
				setOption("BlackBackground", false);
				run("Convert to Mask");
//The ROIs are generated
				run("Analyze Particles...", "size="+conaSize+" pixel circularity="+conaCirc+"-1.00 add include exclude");
//The image is closed
				close();
				open(roi_image);
				roiManager("measure");
				saveAs("Results", inputc+"cona.csv");
				roiManager("Delete");
				run("Clear Results");
				selectWindow("Results");
				run("Close");
				close();
			}
	
		//Now that the nuclear data is collected, Ygg will collected the nuclear indepedent data
		//Since the list of images will be the same, it simply iterates through each image and collects the total ROI pixel data
			if(dapiExtract){
				open(pathway2+"/dna.png");
				run("8-bit");
		// Change this for deviations from default
				setAutoThreshold("Default dark");
		//run("Threshold...");
				setThreshold(dnaMinThreshold, 255);
				setOption("BlackBackground", false);
				run("Convert to Mask");
				run("Watershed");
				run("Analyze Particles...", "size="+nucleusSize+" pixel add include exclude");
				close();
				open(pathway2+"/dna.png");
				roiManager("measure");
				saveAs("Results", inputc+"dna.csv");
				roiManager("Delete");
				run("Clear Results");
				selectWindow("Results");
				run("Close");
				close();
			};
		};
		};
	};
};
selectWindow("ROI Manager");
run("Close");