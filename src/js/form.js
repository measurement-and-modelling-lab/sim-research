const exec = require('child_process').exec;
let spawn = require("child_process").spawn;
const fs = require('fs');
var path = require('path');
var os = require('os');
console.log(os.platform());
var simRunning = false;
var progress = 0;
var numConditions;
var progressBar;
function execute(command, callback) {
    exec(command, (error, stdout, stderr) => {
        callback(stdout);
    });
};
var call = 0;
(function() {
    var run =  document.getElementById("run_button")
    progressBar = document.getElementById("progress_bar")
    run.addEventListener("click", function(){
        // Collect the form data
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
        call = call + 1;
        let testFile = 'test_out'+ call + '.csv';
        var fileToWatch = path.join(__dirname, '../',testFile);
        // build command and run it
        const command = `./simulation ${iterations} ${conditions} ${seeds} >> ${testFile}`
        exec(command, (err, output) => {
            console.log(output);
        });
        // Set this flag so the file will be watched
        simRunning=true;
         //Watch the out file
        fs.watchFile(fileToWatch, function (curr, prev) {
            //console.log(simRunning)
            if(simRunning){
                var i;
                var count = 0;
                //Read current file and get how many lines
                fs.createReadStream(fileToWatch)
                .on('data', function(chunk) {
                    for (i=0; i < chunk.length; ++i)
                    if (chunk[i] == 10) count++;
                })
                .on('end', function() { // After reading update progress bar based on how many lines are there
                    console.log(count);
                    progress = Math.round((((count/2)+1)/numConditions)*100);
                    progressBar.style.width = ""+progress+"%";
                    progressBar.setAttribute('aria-valuenow', progress)
                    progressBar.innerHTML =  ""+progress+"%";
                    console.log(progress);
                });
            }
        });
      });
 })();
