const exec = require('child_process').exec;
const fs = require('fs');
var path = require('path');
var os = require('os');
console.log(os.platform());
var fileToWatch = path.join(__dirname, '../','test_out.csv');
var simRunning = false;
var progress = 0;
var numConditions;
var progressBar;
console.log(fileToWatch);
function execute(command, callback) {
    exec(command, (error, stdout, stderr) => { 
        callback(stdout); 
    });
};

(function() {
    var run =  document.getElementById("run_button")
    progressBar = document.getElementById("progress_bar")
    run.addEventListener("click", function(){
        document.getElementById("progress_bar_div").style.visibility = "visible";
        form = document.getElementById("config_form");
        var conditions = form.conditions.files[0].path;
        var i;
        var count = 0;
        fs.createReadStream(conditions)
            .on('data', function(chunk) {
                for (i=0; i < chunk.length; ++i)
                if (chunk[i] == 10) count++;
            })
            .on('end', function() {
                numConditions = count;
        });
        var iterations = form.iterations.value;
        var seeds = form.seeds.files[0].path;
        console.log(conditions)
        console.log(iterations)
        console.log(seeds)
        // execute('./simulation conditions.csv 12000 seeds.csv', (output) => {
        //     console.log(output);
        // });
        simRunning=true;
      });
 })();
fs.watchFile(fileToWatch, function (curr, prev) {
        console.log(simRunning)
        if(simRunning){
            var i;
            var count = 0;
            fs.createReadStream(fileToWatch)
            .on('data', function(chunk) {
                for (i=0; i < chunk.length; ++i)
                if (chunk[i] == 10) count++;
            })
            .on('end', function() {
                console.log(count);
                progress = Math.round((count/numConditions)*100);
                progressBar.style.width = ""+progress+"%";
                progressBar.setAttribute('aria-valuenow', progress)
                progressBar.innerHTML =  ""+progress+"%";
                console.log(progress);
            });
        }
});