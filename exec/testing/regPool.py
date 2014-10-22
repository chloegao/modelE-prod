import time
import os
import sys
import subprocess
import shlex
import threading
import multiprocessing
import commands
import logging

logger = logging.getLogger('regPool')

#-------------------------------------------------------------------------------
class Worker(multiprocessing.Process):
    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                logger.debug('%s: Exiting' % proc_name)
                self.task_queue.task_done()
                break
            logger.debug('%s: %s' % (proc_name, next_task))
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return

#-------------------------------------------------------------------------------
class Task(object):
    def __init__(self, a):
        self.a = a
    def __call__(self):
        rc = syscmd1(self.a)
        while rc.poll() is None:
            time.sleep(5)
        if rc.returncode !=0:
            logger.debug('%r failed: %s' % (self.a, rc))
        logger.debug('%r is done' % (self.a))
    def __str__(self):
        return '%s' % (self.a)

#-------------------------------------------------------------------------------
class Batch(object):
    def __init__(self, a):
        self.a = a
    def __call__(self):
        rc = sbatchSlurmCmd(self.a)
        while True:
            for job in jobs:
                logger.debug('Monitoring job : ' + job)
                rc=syscmd2('squeue -j '+str(job)+' -t PD,R -h -o %t')
                # If job is done, remove from list
                if rc == '':
                    logger.debug( '...' + job + ' is done')
                    jobs.remove(job)
                else:
                    if 'PD' in rc:
                        logger.debug( '...' + job + ' is pending')
                    else:
                        logger.debug( '...' + job + ' is running')
                    time.sleep(60)
            # If list is empty then we are done
            if not jobs:
                break
    def __str__(self):
        return '%s' % (self.a)

#-------------------------------------------------------------------------------
def runCommands(commands, useBatch):
    # Establish communication queues
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    
    # Start workers
    num_workers = len(commands)
    logger.debug( 'Creating %d workers' % num_workers)
    workers = [ Worker(tasks, results)
                  for i in xrange(num_workers) ]
    for w in workers:
        w.start()
    
    # Enqueue jobs
    num_jobs = len(commands)
    for command in commands:
        if useBatch == 'yes':
            tasks.put(Batch(command))
        else:
            tasks.put(Task(command))

    # Add a poison pill for each worker
    for i in xrange(num_workers):
        tasks.put(None)

    # Wait for all of the tasks to finish
    tasks.join()
  
#-------------------------------------------------------------------------------
#  This function submits a batch job under SLURM and creates a jobs list
def sbatchSlurmCmd(cmd):
    global jobs
    jobs = []
    output = commands.getoutput(cmd % vars())
    # Parse the output from sbatch and append job ID to jobs list
    jobs.append(shlex.split(output)[3])
    return 0

#-------------------------------------------------------------------------------
def syscmd1(cmd):
    return subprocess.Popen(cmd, shell=True)

#-------------------------------------------------------------------------------
def syscmd2(cmd):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    return out




  


