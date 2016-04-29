#!/bin/bash

# A script to run gsaRestApiServer.js with its own Oracle (12.1) client; assumes that Node.js is in the user's $PATH
## source /homes/oracle/ora11setup.sh
## export LD_LIBRARY_PATH=/nfs/production3/ma/home/rpetry/workspace/instantclient:$LD_LIBRARY_PATH
export ATLAS_PROD=$ATLAS_PROD
export PATH=${ATLAS_PROD}/sw/nodejs/node-v0.12.7-linux-x64/bin:$ATLAS_PROD/sw/atlasinstall_prod/R_install/bin/:$PATH

cp ${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/analysis/gsa/api/gsaRestApiServer.js /nfs/production3/ma/home/rpetry/workspace/nodejs/
# gsaRestApiServer.js needs to be run from the same directory in which the modules it requires are installed, c.f. /nfs/production3/ma/home/rpetry/workspace/nodejs/node_modules
node /nfs/production3/ma/home/rpetry/workspace/nodejs/gsaRestApiServer.js > /nfs/production3/ma/home/rpetry/workspace/nodejs/gsaRestApiServer.log 2>&1 &