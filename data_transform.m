clear;
clc;

rootPath = "/Users/haoxiangyang/Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code";
firstsub = ["/ErrorTolerance10","/ErrorTolerance30"];
secondsub = ["/myExperiment30","/myExperiment80"];
lastsub = ["/D1N50","/D1N200","/D1N400","/D1N600","/D1N800","/D1N1000","/D1N1200"];
fileList = ["/Market","/Resource1","/Resource2","/Resource3"];

for i1 = firstsub
    for i2 = secondsub
        for i3 = lastsub
            for i4 = fileList
                try
                    myDataPath = strcat(rootPath,i1,i2,i3,i4,".xlsx");
                    myData = xlsread(myDataPath);
                    csvwrite(strcat(rootPath,i1,i2,i3,i4,".csv"), myData);
                catch
                    a = 1;
                end
            end
        end
    end
end