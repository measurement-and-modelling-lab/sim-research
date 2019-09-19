const exec = require('child_process').exec;
let spawn = require("child_process").spawn;
const fs = require('fs');
var path = require('path');
var os = require('os');
console.log(os.platform());
var simRunning = false;
var progress = 0;
var numConditions;
var numSeeds;
var iterations;
var minSeeds;
var progressBar;
var form;
function execute(command, callback) {
    exec(command, (error, stdout, stderr) => { 
        callback(stdout); 
    });
};

(function() {
    var run_sim = document.getElementById("run_button");
    var gen_seeds = document.getElementById("gen_seeds");
    var conditions_file_element = document.getElementById("conditions_input");
    var seeds_file_element = document.getElementById("seeds_input");
    progressBar = document.getElementById("progress_bar");
    form = document.getElementById("config_form");
    
    
    conditions_file_element.addEventListener("change", function(){
        var conditions = form.conditions.files[0].path;
        fs.createReadStream(conditions){
            .on('data', function(chunk) {
                for (var i=0; i < chunk.length; ++i)
                    if (chunk[i] == 10) numConditions++;
            })
        };
    });
    
    seeds_file_element.addEventListener("change", function(){
        var seeds = form.seeds.files[0].path;
        var numSeeds = 0;
        fs.createReadStream(seeds){
            .on('data', function(chunk) {
                for (i=0; i < chunk.length; ++i)
                    if (chunk[i] == 10) numSeeds++;
            })
        };
    });
    
    gen_seeds.addEventListener("click", function(){
        iterations = form.iterations.value;
        minSeeds = iterations * numConditions;
        let seedsFile = 'seeds(1).csv'
        var fileToWatch = path.join(__dirname, '../',seedsFile);
        // build command and run it 
        const command = `./seed_generator ${minSeeds} >> ${seedFile}`
        
        function run_executable(cmd, callback) {
            exec(cmd, (error, output) => {
                console.log(output);
                callback();
            });
        }
        function finished(){
            console.log("seeds generation finished\n");
        }
        run_executable(command, finished);
    })
        
    run_sim.addEventListener("click", function(){
        // Collect the form data
        document.getElementById("progress_bar_div").style.visibility = "visible";
        
        iterations = form.iterations.value;
        minSeeds = iterations * numConditions;
        if (numSeeds>minSeeds){
            console.log(conditions)
            console.log(iterations)
            console.log(seeds)
            let testFile = 'test_out.csv'
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
                    })
                }
            });
        }else{
            console.log("not enough seeds!\n");
        }
    });
});
