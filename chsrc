##Setup modules for non-interactive PBS jobs
#source /site/modules/default/init/csh 

# User specific environment and startup programs
#module load intel
echo $LOADEDMODULES | grep "PrgEnv-cray" > /dev/null
if ($status == 0) then
  module swap PrgEnv-cray PrgEnv-intel
endif
umask 027
setenv NCARG_ROOT /p/home/wallcraf/pkgs/ncarg/gnu/5.2.1
#set path=($path ~/hycom/ALL/bin ~/bin /site/netcdf/bin ${NCARG_ROOT}/bin .)
set path=($path ~/hycom/ALL/bin ~/bin ${NCARG_ROOT}/bin .)
if ( $?prompt ) then
	set history = 90
	set host    = `uname -n`
	set prompt  = "$host \!> "
#
	set ignoreeof
	set noclobber
	set notify
	setenv PAGER    "less --no-init --hilite-search"
	setenv LESSKEY	$HOME/.less
endif
#
alias cbg  '(csh \!^ |& cat - >& \!^:r.log) &' # csh in bg with output in .log
alias cd   'cd \!*; echo $PWD'                # change directory, and list it
alias cdc  'cd `echo $cwd | awk -f ~wallcraf/bin/cdc.awk`'
                                              # cd to-from $CENTER
alias cdd  'cd `echo $cwd | awk -f ~wallcraf/bin/cdw.awk`'
alias cds  'cd `echo $cwd | awk -f ~wallcraf/bin/cdw.awk`'
alias cdw  'cd `echo $cwd | awk -f ~wallcraf/bin/cdw.awk`'
                                              # cd to-from /work
alias cdp  'cd `echo $cwd $user | awk -f ~wallcraf/bin/cdp.awk`'
                                              # cd to-from projects/hycom
alias colc 'echo "12345678901234567890123456789012345678901234567890123456789012345678901234567890"'
                                              # column count
alias cpu  'cat /proc/cpuinfo'                # type of cpu
alias d2h  '(egrep " ..:..:.. .* 199" \!* /dev/null \\
             | awk -f ~wallcraf/bin/d2h.awk)' # convert dates to wall hrs
alias d2t  '(egrep " ..:..:.. .* 199" \!* /dev/null \\
             | awk -f ~wallcraf/bin/d2t.awk)' # convert dates to wall time
alias dus  '(/usr/bin/du -sk * | sort -n \\
             | awk -f ~wallcraf/bin/dus.awk)' # reformatted du
alias h    'history 20'                       # short history listing
alias idmjdm 'head -n 2 regional.grid.b | sed -e "s/^  *//g" -e "s/ .*//g" | paste -s -d" "'
alias ll   '/bin/ls -laF'                     # list 1 line per file
alias lld  '/bin/ls -lad */.'                 # list directory names
alias mem  'cat /proc/meminfo | head -2'      # how much physical memory
#
alias pg    'less --no-init --hilite-search'    # pg not availabe for Linux
alias pg80  'less --no-init --hilite-search --chop-long-lines'
                                                # pg without line wrap
alias pss  'ps -df \\
            | egrep -v "root |daemon | 00:00:0[0-9] " \\
            | sort'                             # interesting processes
alias psu  'ps -fu ${user} | sort'              # ps of my processes
alias q    '~wallcraf/bin/q_navo'               # qsub with .log file
#alias qa   'qview'                              # PBS batch jobs
#alias qA   'qview'                              # PBS batch jobs
#alias qr   '(qview | grep -v "  [A-Z] ------")' # running PBS batch jobs
#alias qR   '(qview | grep -v "  [A-Z] ------")' # running PBS batch jobs
#alias qu   'qview -u $user'                     # my PBS batch jobs
#alias qU   'qview -u $user'                     # my PBS batch jobs
alias qa   'qstat -a'                         # PBS batch jobs
alias qA   'qstat -a'                         # PBS batch jobs
alias qr   '(qstat -a | grep -v " [A-Z]   --")'   # running PBS batch jobs
alias qR   '(qstat -a | grep -v " [A-Z]   --")'   # running PBS batch jobs
alias qu   'qstat -u $user'                   # my PBS batch jobs
alias qU   'qstat -u $user'                   # my PBS batch jobs
alias pwd  'echo $PWD'                        # pwd hiding /tmp_mnt
alias rm   'rm -ir'                           # always confirm file deletion
alias rsx   'set noglob; \\
             eval `resize -c`; \\
             unset noglob; \\
             set glob; \\
             printenv \\
             | egrep "^COLUM|^LINES"'         # reset TERMCAP for xterm
alias hclog '(egrep "Timer   [12]:|HY_Run" \!* /dev/null \\
              | sed -e "s?......time/call.*??g" -e "s?calls =.*time =?   time =  ?g")'
                                               # HYCOM+CICE wall clock times
alias hlog '(grep "total    calls" \!* /dev/null \\
             | sed -e "s/:.* time =/: /g" -e "s?time/call.*??g" \\
             | awk -f ~wallcraf/bin/hlog.awk)' # print HYCOM wall clock times
alias tlog '(egrep "%[ 	].*pf+" \!* /dev/null \\
             | awk -f ~wallcraf/bin/tlog.awk)' # expand & print "time" lines
alias tlog1 '(egrep "%[ 	].*pf+" \!* /dev/null \\
             | awk -f ~wallcraf/bin/tlog1.awk)' # expand & print "time" lines
alias wtlog '(egrep "^[A-Z].*T 20[0-2][0-9]" \!* /dev/null \\
             | awk -f ~wallcraf/bin/wtlog.awk)' # expand & print "date" lines
alias usage '(/app/bin/show_usage \\
              | awk -f ~wallcraf/bin/usage_erdc.awk)'
                                              # accouting in standard form
alias w80  '~wallcraf/bin/Linux/w80;  \\
            setenv COLUMNS  80; \\
            stty cols  80'                    # put vt100 in  80 column mode
alias w132 '~wallcraf/bin/Linux/w132; \\
            setenv COLUMNS 132; \\
            stty cols 132'                    # put vt100 in 132 column mode

alias QSTAT 'qstat -u grstephe'
alias hgrep 'history | grep'
alias go2d 'cd /p/home/grstephe/hycom/F2Dt025km/expt_02.0'
alias goV 'cd /p/home/vchalama/hycom/F2Dt025km/expt_02.0'
alias go51 'cd /p/home/grstephe/hycom/GLBc0.08/expt_05.1'
alias go61 'cd /p/home/grstephe/hycom/GLBc0.08/expt_06.1'
