pipeline {
    agent any

    environment {
        PATH = "$WORKSPACE/miniconda/bin:$PATH"
    }

    stages {
        stage('setup miniconda') {
            steps {
                sh '''#!/usr/bin/env bash
                wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -nv -O miniconda.sh > /dev/null
                bash miniconda.sh -b -p $WORKSPACE/miniconda > /dev/null
                conda config --set always_yes yes --set changeps1 no > /dev/null
                conda update -q conda  > /dev/null

                # add channels
                conda config --add channels defaults    > /dev/null
                conda config --add channels conda-forge > /dev/null
                conda config --add channels bioconda    > /dev/null

                # create snakemake-workflows env
                conda env create -f envs/conda-bot.yaml > /dev/null
                '''
            }
        }
        stage('Check for updates') {
            steps {
                sh '''#!/usr/bin/env bash
                source $WORKSPACE/miniconda/etc/profile.d/conda.sh
                conda activate miniconda/envs/conda-bot/

                python Jenkins/conda-bot.py $GITHUB_TOKEN
                '''
            }
        }
    }
}
