const exec = require('child_process').exec;
function execute(command, callback) {
    exec(command, (error, stdout, stderr) => { 
        callback(stdout); 
    });
};

(function() {
    var run =  document.getElementById("run_button")
    run.addEventListener("click", function(){
        form = document.getElementById("config_form");
        var conditions = form.conditions.files[0].path;
        var iterations = form.iterations.value;
        var seeds = form.seeds.files[0].path;
        console.log(conditions)
        console.log(iterations)
        console.log(seeds)
        
        // execute('ping -c 4 0.0.0.0', (output) => {
        //     console.log(output);
        // });
      });
 })();