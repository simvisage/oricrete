import os
import subprocess
import errno
import time
import sys

PIPE = subprocess.PIPE

if subprocess.mswindows:
    from win32file import ReadFile, WriteFile
    from win32pipe import PeekNamedPipe
    import msvcrt
else:
    import select
    import fcntl

 
class Popen(subprocess.Popen):
    '''
    This class is a code-snipped from: 
    http://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDgQFjAA&url=http%3A%2F%2Fcode.activestate.com%2Frecipes%2F440554-module-to-allow-asynchronous-subprocess-use-on-win%2Fdownload%2F1%2F&ei=L81SUc6YFdOM4gSWoYHoBQ&usg=AFQjCNEpnjK5yQO8OdMuevJP4gToOp7yvQ&sig2=DgbTCRxkDcR3PSbBOWFdhA&bvm=bv.44342787,d.bGE&cad=rja
    
    ToDo:
    This code can't interpret answers from the cluster,
    so it just waits some time (default 5 sec) and then 
    starting the next input. So it won't catch any errors!!!
    '''
    def recv(self, maxsize = None):
        '''
        Standard output stream from cluster.
        '''
        return self._recv('stdout', maxsize)

    def recv_err(self, maxsize = None):
        '''
        Standard error stream from cluster.
        '''
        return self._recv('stderr', maxsize)

    def send_recv(self, input = '', maxsize = None):
        return self.send(input), self.recv(maxsize), self.recv_err(maxsize)

    def get_conn_maxsize(self, which, maxsize):
        if maxsize is None:
            maxsize = 1024
        elif maxsize < 1:
            maxsize = 1
        return getattr(self, which), maxsize

    def _close(self, which):
        getattr(self, which).close()
        setattr(self, which, None)

    if subprocess.mswindows:
        def send(self, input):
            if not self.stdin:
                return None

            try:
                x = msvcrt.get_osfhandle(self.stdin.fileno())
                (errCode, written) = WriteFile(x, input)
            except ValueError:
                return self._close('stdin')
            except (subprocess.pywintypes.error, Exception), why:
                if why[0] in (109, errno.ESHUTDOWN):
                    return self._close('stdin')
                raise

            return written

        def _recv(self, which, maxsize):
            conn, maxsize = self.get_conn_maxsize(which, maxsize)
            if conn is None:
                return None

            try:
                x = msvcrt.get_osfhandle(conn.fileno())
                (read, nAvail, nMessage) = PeekNamedPipe(x, 0)
                if maxsize < nAvail:
                    nAvail = maxsize
                if nAvail > 0:
                    (errCode, read) = ReadFile(x, nAvail, None)
            except ValueError:
                return self._close(which)
            except (subprocess.pywintypes.error, Exception), why:
                if why[0] in (109, errno.ESHUTDOWN):
                    return self._close(which)
                raise

            if self.universal_newlines:
                read = self._translate_newlines(read)
            return read

    else:
        def send(self, input):
            if not self.stdin:
                return None

            if not select.select([], [self.stdin], [], 0)[1]:
                return 0

            try:
                written = os.write(self.stdin.fileno(), input)
            except OSError, why:
                if why[0] == errno.EPIPE: #broken pipe
                    return self._close('stdin')
                raise

            return written

        def _recv(self, which, maxsize):
            conn, maxsize = self.get_conn_maxsize(which, maxsize)
            if conn is None:
                return None

            flags = fcntl.fcntl(conn, fcntl.F_GETFL)
            if not conn.closed:
                fcntl.fcntl(conn, fcntl.F_SETFL, flags | os.O_NONBLOCK)

            try:
                if not select.select([conn], [], [], 0)[0]:
                    return ''

                r = conn.read(maxsize)
                if not r:
                    return self._close(which)

                if self.universal_newlines:
                    r = self._translate_newlines(r)
                return r
            finally:
                if not conn.closed:
                    fcntl.fcntl(conn, fcntl.F_SETFL, flags)

message = "Other end disconnected!"

def recv_some(p, t = 5, e = 1, tr = 5, stderr = 0):
    '''
    Recieving all answers from cluster in default 5 sec.
    '''
    if tr < 1:
        tr = 1
    x = time.time() + t
    y = []
    r = ''
    pr = p.recv
    if stderr:
        pr = p.recv_err
    while time.time() < x or r:
        r = pr()
        if r is None:
            if e:
                raise Exception(message)
            else:
                break
        elif r:
            y.append(r)
        else:
            time.sleep(max((x - time.time()) / tr, 0))
    return ''.join(y)

def send_all(p, data):
    '''
    Send a row of commands to the shell.
    '''
    while len(data):
        sent = p.send(data)
        if sent is None:
            raise Exception(message)
        data = buffer(data, sent)

def delete_old(p, filename, tail):
    '''
    Delete all old files with the filename and these endings laying on the 
    cluster.
    
    p [Popen]: A active shell.
    filename [string]: Actual file.
    tail [string]: Closing escape sequence (WIN and Unix are different)
    '''
    endings = ['.com', '.log', '.cae', '.dat',
               '.jnl', '.msg', '.odb', '.prt',
               '.sim', '.sta', '.xml', '.stt',
               '.odb_f', '.lck', '.mdl', '.ipm']
    remove = 'rm '
    command_list = []
    for end in endings:
        command_list.append(remove + filename + end + tail)
    for cmd in command_list:
        send_all(p, cmd)

def solve_abaqus(p, filename, tail):
    '''
    Start a abaqus job with the *.inp file.
    
    p [Popen]: A active shell.
    filename [string]: Actual file.
    tail [string]: Closing escape sequence (WIN and Unix are different)
    '''
    send_all(p, 'abaqus job=' + filename + tail)

def open_abaqus(p, tail):
    '''
    Starts the Abaqus CAE.
    
    p [Popen]: A active shell.
    tail [string]: Closing escape sequence (WIN and Unix are different)
    '''
    send_all(p, 'abaqus cae ' + tail)

def connect_cluster(p, login, tail, cluster = 'cluster.rz.rwth-aacheb.de', options = []):
    '''
    Connects to the cluster.
    
    p [Popen]: A active shell.
    login [string]: Working login name.
    tail [string]: Closing escape sequence (WIN and Unix are different)
    cluster [string]: Path of the cluster.
    options [List]: List of option strings.
    '''
    cmd = 'ssh '
    for i in options:
        cmd += i + ' '
    cmd += cluster + ' -l ' + login + tail
    send_all(p, cmd)
    

def upload_file(p, login, path_file, tail, cluster = 'cluster.rz.rwth-aachen.de', path_server = '~/'):
    '''
    Uploads an explicit *.inp file to the cluster.
    
    p [Popen]: A active shell.
    login [string]: Working login name.
    path_file [string]: path of the file, which should be uploaded.
    tail [string]: Closing escape sequence (WIN and Unix are different)
    cluster [string]: Path of the cluster.
    path_server [string]: Path where the file should be uploaded to. Default 
                          is home.
    '''
    cmd = 'scp ' + path_file + ' ' + login + '@' + cluster + ':' + path_server + tail
    send_all(p, cmd)

def download_file(p, login, path_file, tail, path_server, cluster = 'cluster.rz.rwth-aachen.de'):
    '''
    Download an explicit *.inp file from the cluster.
    
    p [Popen]: A active shell.
    login [string]: Working login name.
    path_file [string]: path of the file, where it should be downloaded.
    tail [string]: Closing escape sequence (WIN and Unix are different)
    path_server [string]: Path where the file lays, which should be downloaded to.
    cluster [string]: Path of the cluster.
    '''
    cmd = 'scp ' + login + '@' + cluster + ':' + path_server + ' ' + path_file + tail
    send_all(p, cmd)
    
def close_connection(p, tail):
    '''
    Kills the shell.
    '''
    send_all(p, 'exit' + tail)
    p.kill()
    
def open_shell():
    '''
    Opens a new shell.
    
    Returns a Popen object with the active shell.
    '''
    if sys.platform == 'win32':
        shell, commands, tail = ('cmd', ('dir /w', 'echo HELLO WORLD'), '\r\n')
    else:
        shell, commands, tail = ('sh', ['ssh -Y -t -t cluster-x.rz.rwth-aachen.de -l ms299282', 'echo HELLO WORLD'], '\n')
    a = Popen(shell, stdin = PIPE, stdout = PIPE)
    return a


if __name__ == '__main__':
    tail = '\n'
    #open commandshell 
    a = open_shell()
    

    cluster = 'cluster-x.rz.rwth-aachen.de'
    login = 'ms299282'



    #connect to server
    options = ['-Y',
               '-t',
               '-t']

    connect_cluster(a, login, tail, cluster = cluster, options = options)

    print recv_some(a)

    #remove old datas
    
    #solve new file
    solve_abaqus(a, 'test_name', tail)
    print recv_some(a)

    # close connection
    send_all(a, 'exit' + tail)
    print recv_some(a, e = 0)
    a.wait()

