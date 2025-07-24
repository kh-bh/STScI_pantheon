#!/usr/bin/env python
'''                                                                                                                                                                                            
A wrapper by A. Rest around htc_utils from J. Hunkeler.
http://www.stsci.edu/~jhunk/htcondor/
'''
import os,sys
import subprocess
from htc_utils.bindings import Job, Submit
import argparse 

def makepath(path,raiseError=1):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)


class SubmitCondorBatchClass(Job):
    def __init__(self,filename, ext='job', * args, **kwargs):
       Job.__init__(self,filename, ext=ext,*args, **kwargs)
       self.verbose = 0

       self.cmdlistfilename = self.filename = None

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument('cmdlistfilename',
                            help='filename of list of commands to submit to condor') 
        parser.add_argument('--verbose', '-v', action='count')
        parser.add_argument('--ext',default='job',nargs=1,
                            help='extension of output filename (default=%(default))') 
        parser.add_argument('--machines',default=None,nargs="+",
                            help='list of allowed machines') 
        parser.add_argument('--maxcpus',default=None,nargs="+",type=int,
                            help='How many CPU\'s should be requested?') 
        parser.add_argument('--request_memory',default=None,nargs="+",type=int,
                            help='How much memory should be requested? for DIFFIM (hotpants), ask for 8192 (default=%(default))') 
        parser.add_argument('--jobfilename',default='auto',nargs=1,
                            help='specify jobfilename which is submitted to condor. If \'auto\', then it is input filename with the extension ext (default=%(default)s)') 
        parser.add_argument('--skiplog',help="Don't save the log files for the cluster of jobs. Otherwise the log files are saved into <basedir>/logs/<basename> directory, with the <basedir>/<basename>=jobfilename) ",action="store_true",default=False)
        parser.add_argument('--logdir',default=None,type=str,
                            help='specify the logdir. if not specified, then the logdir is <basedir>/logs/<basename> directory, with the <basedir>/<basename>=jobfilename. Note that you can skip the log files with --skiplog. (default=%(default)s)')
        parser.add_argument('--stdout2log',help="Save the stdout of each job to its own log file (default=%(default))",action="store_true",default=False)
        parser.add_argument('--stderr2log',help="Save the stderr of each job to its own log file (default=%(default))",action="store_true",default=False)
        return(parser)

    def setoptions(self,args):
        ''' set the attributes of this job from the command line args, which are parsed with parse_args()
        '''

        self.verbose = args.verbose

        self.cmdlistfilename =  os.path.abspath(args.cmdlistfilename)
        
        # set the filename
        if args.jobfilename.lower()=='auto':
            self.filename = '.'.join([self.cmdlistfilename, args.ext])
            if os.path.isfile(self.filename) and os.path.samefile(self.cmdlistfilename,self.filename):
                raise RuntimeError('job filename is the same than inputfilename, refusing to overwrite input filename %s!' % (self.cmdlistfilename))
        else:
            self.filename = os.path.abspath(args.jobfilename)
        
        if not args.skiplog:
            print(f'CCCCCCC {args.logdir}')
            if args.logdir is None:
                (basedir,basename)=os.path.split(self.filename)
                logdir = f'{basedir}/logs/{basename}'
            else:
                logdir = args.logdir 
            print(f'LOG DIRECTORY {logdir}')
            makepath(logdir)
            self.logging(logdir, create=True)    
        else:
            print('skipping logging...')

        if args.machines is not None:
            if isinstance(args.machines,str):
                machines = [args.machines]
            elif isinstance(args.machines,list):
                machines = args.machines
            else:
                raise RuntimeError(f'machines should be string or list, {args.machines} given!')

            counter=0
            s=''
            if len(machines)>1:s='(';
            for machine in machines:
                if counter>0: s+=' || '
                s+=f'machine == "{machine}"'
                counter+=1
            if len(machines)>1:s+=')';
            print('Adding requirement:',s)
            self.attr('requirements', s) 

        # if args.maxcpus is not None:
        # s = f'request_cpus == {args.maxcpus}'
        # print('Adding requirement:',s)
        # self.attr('requirements', s)
        if args.maxcpus is not None:      
            self.attr('request_cpus',args.maxcpus)

        if args.request_memory is not None:      
            self.attr('request_memory',args.request_memory)



        return(0)

    def parse_cmds(self,cmdlist):
        fullexe0 = None
        exe0=None

        counter=0

        # loop through the cmd lines
        for c in cmdlist:
            c = c.strip()
            (exe,arguments)=c.split(None,1)

            # make sure all executable are the same
            if exe0 == None:
                exe0=exe

                # ok, not the best, but compatible with python 2.7 and 3: use which to get full executable name
                my_env = os.environ.copy()
                out = subprocess.Popen('which %s' % exe, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True, shell=True)
                cmd_in, cmd_out = out.stdin, out.stdout

                output = cmd_out.readlines()
                if len(output)<1:
                    raise RuntimeError('Could not find the executable %s' % exe)
                fullexe0 =output[0].strip().decode('UTF-8')
                
                # just one more test, shoudl be changed to check if it is an executable
                if not os.path.isfile(fullexe0):
                    raise RuntimeError('executable %s does not exist!' % fullexe0)
 
                self.attr('executable', fullexe0)
            else:
                if exe0 != exe:
                    raise RuntimeError('Can only run a list of commands with the same executable %s, but %s is alo passed as executable!' % (exe0,exe))
                # all good!!

            # add the cmd to the queue!
            if not SubmitCondorBatch.verbose is None and self.verbose>1: print('job %d: Adding command arguments to queue: %s ' % (counter,arguments))
            self.subattr('arguments', arguments)
            self.subattr('queue')

            counter +=1
        if self.verbose: print('%d jobs saved into %s' % (counter,self.filename))
        return(0)

    def parse_cmdlistfile(self,cmdlistfilename):
        cmdlist = open(cmdlistfilename,'r').readlines()
        return(self.parse_cmds(cmdlist))

    def submitjobs2condor(self):
        # Pass our job to the Submit class
        self.submit = Submit(self)

        # Send our job to the cluster
        self.submit.execute()



if __name__ == '__main__':
    SubmitCondorBatch=SubmitCondorBatchClass('NONE')

    # parse the arguments
    parser = SubmitCondorBatch.add_options()
    args = parser.parse_args()

    # go through the command line options
    SubmitCondorBatch.setoptions(args)

    if SubmitCondorBatch.verbose is not None and SubmitCondorBatch.verbose>0:
        print('cmd list filename:',SubmitCondorBatch.cmdlistfilename)
        print('Jobfilename:',SubmitCondorBatch.filename)

    # parse the input file
    SubmitCondorBatch.parse_cmdlistfile(SubmitCondorBatch.cmdlistfilename)

    # save the condor job file
    if SubmitCondorBatch.verbose is not None and SubmitCondorBatch.verbose: print('Saving',SubmitCondorBatch.filename)
    SubmitCondorBatch.commit()

    #submit the file
    SubmitCondorBatch.submitjobs2condor()

