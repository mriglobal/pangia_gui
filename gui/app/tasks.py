#import json
import sys
import os
import shutil
import time
import datetime
#import pysam
import ast
#from subprocess import Popen, PIPE
from pathlib import Path
#from flask import render_template
from rq import get_current_job
from app import create_app, db
from app.models import User, Task, Results, ResultsToItem, Item, Category
from werkzeug.utils import secure_filename

app = create_app()
app.app_context().push()

basedir = os.path.abspath(os.path.dirname(__file__))

def _set_task_progress(progress):
    job = get_current_job()
    if job:
        job.meta['progress'] = progress
        job.save_meta()
        task = Task.query.get(job.get_id())
        # task.user.add_notification('task_progress', {'task_id': job.get_id(), 'progress': progress})
        if progress >= 100:
            task.complete = True
        db.session.commit()


def _set_real_time_results(new_results_id):
    job = get_current_job()
    if job:
        task = Task.query.get(job.get_id())
        opts = ast.literal_eval(task.options)
        res_dict = {'results_id': new_results_id}
        new_options = {**opts, **res_dict}
        task.options = '{}'.format(new_options)
        db.session.commit()


def remove_dir(directory):
    directory = Path(directory)
    for item in directory.iterdir():
        if item.is_dir():
            remove_dir(item)
        else:
            item.unlink()
    directory.rmdir()


def run_pangia(user_id, args):
    try:

        _set_task_progress(0)

        # Make directories
        os.makedirs(args['tmp_dir'], exist_ok=True)

        # Run preprocessing
        if args['preprocess'] == 'y':
            fastq_split = args['fastq_loc'].split(' ')
            # Paired fastqs
            if len(fastq_split) == 2:
                fastq1_processed = '{}/{}.fastq'.format(args['tmp_dir'], Path(fastq_split[0]).stem)
                fastq2_processed = '{}/{}.fastq'.format(args['tmp_dir'], Path(fastq_split[1]).stem)
                fastp_input = '-i {} -I {}'.format(fastq_split[0], fastq_split[1])
                fastp_output = '-o {} -O {}'.format(fastq1_processed, fastq2_processed)
                pangia_input = '{} {}'.format(fastq1_processed, fastq2_processed)
            # Single fastq
            else:
                fastq1_processed = '{}/{}.fastq'.format(args['tmp_dir'], Path(fastq_split[0]).stem)
                fastp_input = '-i {}'.format(fastq_split[0])
                fastp_output = '-o {}'.format(fastq1_processed)
                pangia_input = fastq1_processed

            if args['trim_polya'] == 'y':
                fastp_trim_polya = ''
            else:
                fastp_trim_polya = ' -A '
            fastp_cmd = 'fastp {}{} -q {} -e {} -l {} -n {} -Y {} -f {} -t {} -h {} -j {} {}'.format(
                fastp_input,
                fastp_trim_polya,
                args['trim_quality'],
                args['avg_qual_cutoff'],
                args['min_read_len'],
                args['n_base_cutoff'],
                args['low_comp_filter'],
                args['cut_5_prime'],
                args['cut_3_prime'],
                '{}/fastp.html'.format(args['tmp_dir']),
                '{}/fastp.json'.format(args['tmp_dir']),
                fastp_output
            )
            # Run Fastp
            fastp_results = os.system(fastp_cmd)

        # No preprocessing
        else:
            pangia_input = args['fastq_loc']

        _set_task_progress(20)

        # Run PanGIA
        # -asl seed length (int) default 40
        # -ams alignment score (int) default 60
        # -sm score method (string) default standalone
        # -ms min score (float) default 0
        # -mr min_read_count (int) default 10
        # -mb min_read_rsnb (int) default 2.5?
        # -ml min_linear_len (int) default 200
        # -mc min_genome_cov (float) default 0.004
        # -md min_depth (float) default 0.01
        # -mrd min_rs_depth (float) default 0.0009
        # -pd pathogen_discovery (y/n) default n
        # --annoy run_annoy (y/n) default n
        # --markers run_tmark (y/n) default n

        pangia_pd = ''
        pangia_annoy = ''
        pangia_tmark = ''
        if args['pathogen_discovery'] == 'y':
            pangia_pd = ' -pd'
        if args['run_annoy'] == 'y':
            pangia_annoy = ' --annoy'
        if args['run_tmark'] == 'y':
            pangia_tmark = ' --markers'

        pangia_cmd = 'python {}pangia.py -i {} -db {}*.fa -st {} -t {} {}{}{} -asl {} -ams {} -sm {} -ms {} -mr {} -mb {} -ml {} -mc {} -md {} -mrd {} -o {}'.format(
            args['pangia_loc'],
            pangia_input,
            args['pangia_db'],
            args['seq_type'],
            args['threads'],
            pangia_pd,
            pangia_annoy,
            pangia_tmark,
            args['seed_length'],
            args['min_align_score'],
            args['score_method'],
            args['min_score'],
            args['min_read_count'],
            args['min_read_rsnb'],
            args['min_linear_len'],
            args['min_genome_cov'],
            args['min_depth'],
            args['min_rs_depth'],
            args['tmp_dir']
        )
        results = os.system(pangia_cmd)
        _set_task_progress(80)

        # Save Results <-------- Saving PanGIA parameters as a string dict. Think of a way to move to a sql db so it's queryable
        new_results = Results(name=args['runname'], description=args['description'], results_type='pangia',
                              user_id=user_id, results_param=str(args))
        db.session.add(new_results)
        db.session.flush()

        results_folder = args['tmp_dir'].replace('/tmp', '/{}'.format(new_results.id))
        new_results.path = results_folder

        for i in args['items_run']:
            new_rtoi = ResultsToItem(results_id=new_results.id, item_id=i)
            db.session.add(new_rtoi)
        db.session.commit()
        _set_task_progress(90)

        # Move Files
        os.makedirs(results_folder, exist_ok=True)
        os.rename('{}/{}.report.tsv'.format(args['tmp_dir'], args['file_stem']), '{}/{}.report.tsv'.format(results_folder, args['file_stem']))
        os.rename('{}/{}.pangia.log'.format(args['tmp_dir'], args['file_stem']), '{}/{}.pangia.log'.format(results_folder, args['file_stem']))
        if args['preprocess'] == 'y':
            os.rename('{}/fastp.html'.format(args['tmp_dir']), '{}/fastp.html'.format(results_folder))
            os.rename('{}/fastp.json'.format(args['tmp_dir']), '{}/fastp.json'.format(results_folder))

        # Cleanup Files
        shutil.rmtree(args['tmp_dir'])

        _set_task_progress(100)

    except:
        _set_task_progress(100)
        app.logger.error('PanGIA Error', exc_info=sys.exc_info())


def run_real_time(user_id, args):
    try:
        _set_task_progress(10)

        # Make directory if it doesn't exist
        if not os.path.exists(args['fastq_dir']):
            os.makedirs(args['fastq_dir'])

        have_been_run = []

        # Infinate Loop through fastqs in this folder
        while True:

            fastqs = [x for x in os.listdir(args['fastq_dir']) if Path(x).suffix == '.fastq' and x not in have_been_run]
            if len(fastqs) == 0:
                print('No more fastq files to run')
                time.sleep(20)
            else:
                # Loop through all of the fastq files
                for fastq in fastqs:
                    fastq_path = '{}/{}'.format(args['fastq_dir'], Path(fastq).name)
                    fastq_name = Path(fastq).name

                    # Loop through all of the databases and run minimap
                    mmis = [x for x in os.listdir(args['pangia_db']) if Path(x).suffix == '.mmi']
                    for mmi in mmis:

                        # MINIMAP
                        this_sam = '{}/{}{}.sam'.format(
                            Path(args['fastq_dir']).resolve(),
                            Path(fastq).stem,
                            Path(mmi).name
                        )
                        minimap_cmd = 'minimap2 -a --sam-hit-only -t {} {}/{} {} -o {}'.format(
                            args['threads'],
                            Path(args['pangia_db']),
                            Path(mmi).name,
                            fastq_path,
                            this_sam
                        )
                        #print(minimap_cmd)
                        minimap = os.system(minimap_cmd)

                        # SAMTOOLS
                        this_bam = '{}/{}{}.bam'.format(
                            Path(args['fastq_dir']).resolve(),
                            Path(fastq).stem,
                            Path(mmi).name
                        )
                        samtools_cmd = 'samtools view -bh -T {}{} {} > {}'.format(
                            args['pangia_db'],
                            Path(mmi).stem,
                            this_sam,
                            this_bam
                        )
                        #print(samtools_cmd)
                        samtools = os.system(samtools_cmd)

                    # Remove small bams files (These were failing to merge)
                    for bam in os.listdir(args['fastq_dir']):
                        if Path(bam).suffix == '.bam':
                            if Path('{}/{}'.format(args['fastq_dir'], bam)).stat().st_size < 1000:  # Less than 1kb
                                remove = os.system('rm {}/{}'.format(args['fastq_dir'], bam))

                    # Merge bam files and convert back to sam file
                    merge_cmd = 'samtools merge -uf {}/master.bam {}/*.bam'.format(args['fastq_dir'], args['fastq_dir'])
                    #print(merge_cmd)
                    merge = os.system(merge_cmd)
                    bam_to_sam_cmd = 'samtools view {}/master.bam > {}/master.sam'.format(args['fastq_dir'], args['fastq_dir'])
                    #print(bam_to_sam_cmd)
                    bam_to_sam = os.system(bam_to_sam_cmd)
                    gawk = os.system('bash {}/rt_gawk.sh {}'.format(basedir, args['fastq_dir']))  # Use gawk.sh outputs master.gawk.sam that will be run through pangia
                    mv_master = os.system('mv -f {}/master.bam {}/master.bak'.format(args['fastq_dir'], args['fastq_dir']))  # Move master bam so it doesn't get deleted

                    # Check for decision tree files
                    dt_db = '{}/PANMLmodels'.format(args['pangia_db'])
                    if os.path.exists(dt_db):
                        run_dt = ' -panml'
                    else:
                        run_dt = ''
                        print('No decision tree databases found.')

                    # Run PanGIA
                    pangia_cmd = 'python {}pangia.py{} -st nanopore -da -db {}*.fa -s {}/master.gawk.sam -o {}/master.results -t {}'.format(
                        args['pangia_loc'],
                        run_dt,
                        args['pangia_db'],
                        args['fastq_dir'],
                        args['fastq_dir'],
                        args['threads']
                    )
                    #print(pangia_cmd)
                    pangia_run = os.system(pangia_cmd)

                    # Remove all bam and sam files
                    remove = os.system('rm -f {}/*.bam'.format(args['fastq_dir']))
                    remove = os.system('rm -f {}/*.sam'.format(args['fastq_dir']))
                    mv_master = os.system('mv -f {}/master.bak {}/last_run.bam'.format(args['fastq_dir'], args['fastq_dir']))

                    # Add to list of finished fastqs
                    have_been_run.append(fastq)

                    # Add Results to DB
                    if len(have_been_run) == 1:
                        new_results = Results(name=args['runname'], description=args['description'], results_type='real_time',
                                              user_id=user_id, results_param=str(args))
                        db.session.add(new_results)
                        db.session.flush()
                        print('Result ID: {}'.format(new_results.id))

                    # Combine fastqs and move them
                    project_folder = '{}/{}'.format(args['upload_dir'], args['project'])
                    results_folder = '{}/results/{}'.format(project_folder, new_results.id)
                    new_fq_name = '{}/{}/{}.fastq'.format(args['upload_dir'], args['project'],  secure_filename(args['runname']))
                    if not os.path.exists(results_folder):
                        create_folder = os.makedirs(results_folder)
                    combine_fq = os.system('cat {}/*.fastq > {}'.format(args['fastq_dir'], new_fq_name))

                    # Move results
                    mv_results = os.system('cp -f {}/master.results/master.gawk.report.tsv {}/{}.report.tsv'.format(args['fastq_dir'], results_folder, secure_filename(args['runname'])))

                    # TIMESTAMPED RESULTS
                    now = datetime.datetime.now()
                    time_results = os.system('cp -f {}/master.results/master.gawk.report.tsv {}/{}.{}.report.tsv'.format(args['fastq_dir'], args['fastq_dir'], secure_filename(args['runname']), now.strftime('%Y_%m_%d_%H:%M:%S')))

                    mv_log = os.system('cp -f {}/master.results/master.gawk.pangia.log {}/{}.pangia.log'.format(args['fastq_dir'], results_folder, secure_filename(args['runname'])))

                    # Add Item to DB
                    if len(have_been_run) == 1:
                        project = Category.query.filter_by(slug=args['project']).first()
                        new_item = Item(name='{}.fastq'.format(secure_filename(args['runname'])),
                            item_path=new_fq_name, category_id=project.id, meta_template_id=0)
                        db.session.add(new_item)
                        db.session.flush()
                        print('Item ID: {}'.format(new_item.id))

                        new_results.path = '{}/results/{}/'.format(project_folder, new_results.id)
                        new_rtoi = ResultsToItem(results_id=new_results.id, item_id=new_item.id)
                        db.session.add(new_rtoi)
                        db.session.commit()

                        # Add results id to the task
                        _set_real_time_results(new_results.id)

    except:
        _set_task_progress(100)
        app.logger.error('PanGIA Error', exc_info=sys.exc_info())


def real_time_finish(user_id, args):
    # Cat all fastqs
    merge_fq = os.system('cat {}/*.fastq > final.fastq'.format(args['fastq_dir']))

    # Add fastq to project folder

    # Add item to database

    # Run one last PanGIA analysis on the final fastq