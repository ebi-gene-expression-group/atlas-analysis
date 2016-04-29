// Based on https://github.com/sbalagop/neo/blob/master/nserver.js
var express = require('express');
var app = express();
var morgan = require('morgan');
var fs = require('fs');
var lineByLineReader = require('line-by-line');
var hashMap = require('hashmap');
var bodyParser = require('body-parser');
var FileStreamRotator = require('file-stream-rotator');
var csv2json = require('csvtojson').Converter;
var childProcess = require('child_process').execSync;
var validator = require('express-validator');
// Cluster implementation adopted from: http://www.sitepoint.com/how-to-create-a-node-js-cluster-for-speeding-up-your-apps/
var cluster = require('cluster');

// Use body parser to parse JSON body
app.use(bodyParser.json());

// Allow access to the X-Forwarded-* header fields - to recover IP address of the client behind the EBI load balancer
// c.f. http://stackoverflow.com/questions/7139369/remote-ip-address-with-node-js-behind-amazon-elb 
app.enable('trust proxy');

// Additional custom validator for checking if the number of genes is less tha the maximum allowed
app.use(validator({
 customValidators: {
     numGenesNoMoreThan: function(value, max) {
         return value.split(" ").length <= max;
    },
 }
}));

var workingDir="/tmp";
var dataDir = "/ebi/ftp/pub/databases/arrayexpress/data/atlas/gsa";
// Based on https://github.com/expressjs/morgan (log file rotation example)
var logDir = "$ATLAS_PROD/atlas_gsa_api/logs";
// ensure log directory exists
fs.existsSync(logDir) || fs.mkdirSync(logDir)
// create a rotating write stream
    var accessLogStream = FileStreamRotator.getStream({
	    filename: logDir + '/access-%DATE%.log',
	    frequency: 'daily',
	    date_format: 'YYYYMMDD',
	    verbose: true
	});
// setup the logger
app.use(morgan('combined', {stream: accessLogStream}))

// Global variables
// The map containing expAcc->contrastID->contrastTitle mapping
var exp2ContrastId2Title = new hashMap();
// The file containing expAcc->contrastID->contrastTitle mapping that will be loaded into exp2ContrastId2Title
var lr = new lineByLineReader(dataDir+'/contrastTitles.tsv');

// Default list of allowed characters for string parameters, e.g. ORGANISM or EFO_TERM
var defaultWhiteList = '\\(\\)\\_\\.\\+\\-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ ';

// Function definitions
// Log any errors
function logError(err) {
    if (err) {
        console.error(err.message);
    }
}

// Return data in JSON format: format (json or tsv);
function returnResults(res, err, resultS, format, errMessage) { 
    if (format === "json") { 
	res.set('Content-Type', 'application/json');
	if (err) {
	    res.status(500).send(JSON.stringify({ status: 500, message: errMessage }));
	} else {
	    var converter = new csv2json({delimiter:"\t"}); 
	    //end_parsed will be emitted once parsing finished
	    converter.on("end_parsed", function (jsonArray) {
		    // console.log(jsonArray); // debug of the resulting jsonarray
		});
	    converter.fromString(resultS, function(err,result){
		    res.status(200).send(JSON.stringify(result));
		});
	}
    } else if (format === "tsv") {
        res.set('Content-Type', 'text/plain');
	if (err) {
            res.status(500).send(errMessage);
        } else {
	    res.status(200).send(resultS);
	}
    } else {
	res.set('Content-Type', 'text/plain');
	res.status(500).send("Unrecognised data format:" + format);
	logError("Unrecognised format:" + format);
    }
} 

// Validate all possible parameters
function paramsValid(req, res) {
    if (req.params.EXPACC)
	req.checkParams("EXPACC","Invalid experiment accession").isWhitelisted(defaultWhiteList);
    if (req.params.CONTRASTID)
        req.checkParams("CONTRASTID","Invalid comparison identifier").isWhitelisted(defaultWhiteList);
    if (req.params.FORMAT)
	req.checkParams("FORMAT","The only formats allowed are tsv and json").isIn(['tsv','json']);
    if (req.params.ORGANISM)
	req.checkParams("ORGANISM","Invalid organism name").isWhitelisted(defaultWhiteList);
    if (req.params.GENE_IDS) {
	req.checkParams("GENE_IDS","Invalid gene identifiers").isWhitelisted(defaultWhiteList);
	req.checkParams('GENE_IDS', 'The number of gene identifiers must be no more than 100').numGenesNoMoreThan(100);
    }

    errors=req.validationErrors();
    if (errors) {
	res.send(errors);
	return false;
    } else {
	return true;
    }
}

// ********* REST API calls start here ***********  

// e.g. http://plantain:3001/getOrganisms 
app.get('/getOrganisms', function (req, res) {
        "use strict";
        if (!paramsValid(req, res)) return;

	try {
	    var lsOrganisms = childProcess("ls " + dataDir + "/*.po | sed 's|\.po$||g' | awk -F\"/\" '{print $NF}'");
            res.set('Content-Type', 'text/plain');
	    res.status(200).send(lsOrganisms);
	} catch (ex) {
	    returnResults(res, ex, null, req.params.FORMAT, "Error retrieving available organisms");
	    return;
	}
    });



// e.g. http://plantain:3001/getContrastTitle/E-TABM-90/g4_g3
app.get('/getContrastTitle/:EXPACC/:CONTRASTID', function (req, res) {
	"use strict";

	if (!paramsValid(req, res)) return;

	var result;
	var err;
	if (exp2ContrastId2Title.has(req.params.EXPACC))
	    result = exp2ContrastId2Title.get(req.params.EXPACC).get(req.params.CONTRASTID);
	res.set('Content-Type', 'text/plain');
	res.status(200).send(result);
    });

// Retrieve overlapping comparisons with :GENE_IDS in :ORGANISM
// e.g. http://plantain:3001/tsv/getOverlappingComparisons/arabidopsis_thaliana/AT1G48030 AT1G53240 AT2G17130 AT2G20420 AT2G44350 AT2G47510 AT3G09810 AT3G15020 AT3G17240 AT3G27380 AT3G55410 AT3G60100 AT4G26910 AT4G35260 AT4G35650 AT4G35830 AT5G03290 AT5G08300 AT5G23250 AT5G40650 AT5G50950 AT5G55070 AT5G65165 AT5G65750 AT5G66760
app.get('/:FORMAT/getOverlappingComparisons/:ORGANISM/:GENE_IDS', function (req, res) {
	"use strict";

	if (!paramsValid(req, res)) return;

        var hrTime = process.hrtime();
        var timeInMicros = hrTime[0] * 1000000 + hrTime[1] / 1000;
        var resultFile = workingDir+"/gsa."+req.params.ORGANISM+timeInMicros;
	try {
	    var gsaCmd="gsa_run.R  --db " + dataDir + "/" + req.params.ORGANISM + ".po --gs '"+ req.params.GENE_IDS +"'  --pvalue 0.05 --out " + resultFile + " -c 4";
	    console.log("Running command: " + gsaCmd);
	    var start = new Date().getTime() / 1000;
	    var gsaRun = childProcess(gsaCmd);
	    var duration = Math.floor((new Date().getTime() / 1000) - start);
	    console.log("Completed successfully in " + duration + " seconds command: " + gsaCmd);
	    // process.stdout.write(gsaRun); // output from gsa_run.R of the above command
	} catch (ex) {
	    returnResults(res, ex, null, req.params.FORMAT, "Error retrieving overlapping Atlas comparisons for organism: '" + req.params.ORGANISM + "' and genes: '" + req.params.GENE_IDS + "'");
	    return;
	}
	// Convert resultFile to the final version to be output to the user and put it in outArr
        var outArr = [];
	var lr = new lineByLineReader(resultFile + ".tsv");
	lr.on('error', function (err) {
		return console.log(err);
	    });

	lr.on('line', function (line) {
	    var arr = line.split("\t");
	    var exp = arr[0];
	    if (exp == "exp") { // Replace the first line
		if (req.params.FORMAT == "tsv" ) {
		    outArr.push("# Enrichment (Fisher-exact, FDR=0.01) of the following gene set across differentially expressed genes in each " + req.params.ORGANISM + " comparison in Expression Atlas: ");
		    outArr.push("# '" + req.params.GENE_IDS + "'");
		}
		outArr.push("EXPERIMENT\tCOMPARISON_ID\tP-VALUE\tOBSERVED\tEXPECTED\tADJUSTED P-VALUE\tEFFECT SIZE\tCOMPARISON_TITLE\tEXPERIMENT_URL");
	    } else {
		    var arrExp = exp.split(":");
		    var expAcc = arrExp[0].split("_")[0];
		    var contrastId = arrExp[3];
		    var expURL = "http://www.ebi.ac.uk/gxa/experiments/"+expAcc+"?queryFactorValues="+contrastId+"&_specific=on";
		    var outLineArr = [];
		    outLineArr.push(expAcc); // expAcc
		    outLineArr.push(contrastId); // contrastId
		    outLineArr.push(arr[1]); // pval
		    outLineArr.push(arr[2]); // observed
		    outLineArr.push(arr[3]); // expected
		    outLineArr.push(arr[4]); // adjPvalue
		    outLineArr.push(arr[5]); // effectSize
                    if (exp2ContrastId2Title.has(expAcc)) {
                        outLineArr.push(exp2ContrastId2Title.get(expAcc).get(contrastId)); // contrastTitle
                    } 
		    outLineArr.push(expURL); // experiment Atlas URL
		    outArr.push(outLineArr.join("\t"));
		}
	    });
	lr.on('end', function (err) {
		// All lines are read, file is closed now - output outArr as the final result in format: req.params.FORMAT
		returnResults(res, err, outArr.join("\n"), req.params.FORMAT, null); 

		// Delete temporary files.
		fs.unlinkSync(resultFile+".tsv");
	    });
    });

if(cluster.isMaster) {
    var numWorkers = require('os').cpus().length;

    console.log('Master cluster setting up ' + numWorkers + ' workers...');

    for(var i = 0; i < numWorkers; i++) {
	cluster.fork();
    }

    cluster.on('online', function(worker) {
	    console.log('Worker ' + worker.process.pid + ' is online');
	});

    cluster.on('exit', function(worker, code, signal) {
	    console.log('Worker ' + worker.process.pid + ' died with code: ' + code + ', and signal: ' + signal);
	    console.log('Starting a new worker');
	    cluster.fork();
	});
} else {
    app.all('/*', function(req, res) {res.send('process ' + process.pid + ' says hello!').end();});

    var server = app.listen(3001, function () {
	"use strict";

	var host = server.address().address,
	port = server.address().port;

	// Read expAcc->contrastID->contrastTitle mapping into exp2ContrastId2Title
	lr.on('error', function (err) {
		return console.log(err);
	    });
	lr.on('line', function (line) {
		var arr = line.split("\t");
		var expAcc=arr[0];
		var contrastID=arr[1];
		var contrastTitle=arr[2];
		if (!exp2ContrastId2Title.has(expAcc)) {
		    exp2ContrastId2Title.set(expAcc, new hashMap());
		}
		exp2ContrastId2Title.get(expAcc).set(contrastID,contrastTitle);
	    });

	console.log('Read in expAcc->contrastID->contrastTitle mapping successfully'); 

	console.log(' Server is listening at http://%s:%s', host, port);
	});
}
	