run("8-bit");
run("Set Scale...", "distance=0.7881 known=10 unit=µm");
setAutoThreshold("Default");
//run("Threshold...");
setThreshold(0, 105);
//setThreshold(0, 105);
setOption("BlackBackground", false);
run("Make Binary", "thresholded remaining black");
run("Analyze Particles...", "size=10-60 circularity=0.40-1.00 show=Outlines summarize add");
