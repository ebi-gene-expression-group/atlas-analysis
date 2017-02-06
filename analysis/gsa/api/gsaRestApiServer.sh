#!/bin/bash

# A script to run gsaRestApiServer.js
export ATLAS_PROD=$ATLAS_PROD
export PATH=${ATLAS_PROD}/sw/nodejs/node-v0.12.7-linux-x64/bin:$ATLAS_PROD/sw/atlasinstall_prod/R_install_3.3.0_rh7/bin/:$PATH

cp ${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/analysis/gsa/api/gsaRestApiServer.js $ATLAS_PROD/atlas_gsa_api
# gsaRestApiServer.js needs to be run from the same directory in which the modules it requires are installed, c.f. /nfs/production3/ma/home/rpetry/workspace/nodejs/node_modules
node $ATLAS_PROD/atlas_gsa_api/gsaRestApiServer.js > $ATLAS_PROD/atlas_gsa_api/gsaRestApiServer.log 2>&1 &
