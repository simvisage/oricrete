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
    def recv(self, maxsize = None):
        return self._recv('stdout', maxsize)
    
    def recv_err(self, maxsize = None):
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
    while len(data):
        sent = p.send(data)
        if sent is None:
            raise Exception(message)
        data = buffer(data, sent)
        
def delete_old(p, filename, tail):
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
    send_all(p, 'abaqus job=' + filename + tail)
    
def load_abaqus(p, filename, tail):
    send_all(p, 'abaqus cae ' + tail)

def connect_cluster(p, login, tail, cluster = 'cluster.rz.rwth-aacheb.de', options = []):
    cmd = 'ssh '
    for i in options:
        cmd += i + ' '
    cmd += cluster + ' -l ' + login + tail
    send_all(p, cmd)

def upload_file(p, login, path_file, tail, cluster = 'cluster.rz.rwth-aachen.de', path_server = '~/'):
    cmd = 'scp ' + path_file + ' ' + login + '@' + cluster + ':' + path_server + tail
    send_all(p, cmd)
    
def download_file(p, login, path_file, tail, path_server, cluster = 'cluster.rz.rwth-aachen.de'):
    cmd = 'scp ' + login + '@' + cluster + ':' + path_server + ' ' + path_file + tail
    send_all(p, cmd)
    

if __name__ == '__main__':
    if sys.platform == 'win32':
        shell, commands, tail = ('cmd', ('dir /w', 'echo HELLO WORLD'), '\r\n')
    else:
        shell, commands, tail = ('sh', ['ssh -Y -t -t cluster-x.rz.rwth-aachen.de -l ms299282', 'echo HELLO WORLD'], '\n')
    
    
    #open commandshell 
    a = Popen(shell, stdin = PIPE, stdout = PIPE)
    
    
    cluster = 'cluster-x.rz.rwth-aachen.de'
    login = 'ms299282'
    
#    upload_file(a, login, '/home/matthias/beamtest.py', tail)
    download_file(a, login, '/home/matthias/', tail, '~/beamtest.py')
    print recv_some(a)
    #connect to server
    options = ['-Y',
               '-t',
               '-t']
    
    connect_cluster(a, login, tail, cluster = cluster, options = options)
    

    
    print recv_some(a)
    
    #remove old datas
    delete_old(a, 'beam', tail)
    
    #solve new file
    solve_abaqus(a, 'beam', tail)
    print recv_some(a)
    
    # close connection
    send_all(a, 'exit' + tail)
    print recv_some(a, e = 0)
    a.wait()

