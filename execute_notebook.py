"""Run the standalone notebook in a NEW real IPython/Jupyter kernel via nbclient.

The working directory is empty and outside the repository. The default is an
external kernel process; --kernel inprocess uses the official ipykernel local
implementation where socket binding is restricted. Both use Jupyter messages. Only the
embedded notebook sources are used. Exceptions stop execution; saved notebook
outputs and a separate execution record identify what was actually run.
"""
import argparse
import hashlib
import os
import json
from pathlib import Path
import platform
import sys
import tempfile
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]


def execute(kernel_mode="process"):
    import nbformat
    from nbclient import NotebookClient
    from jupyter_client.manager import KernelManager
    from jupyter_client.kernelspec import KernelSpec,KernelSpecManager
    class CurrentPython(KernelSpecManager):
        def get_kernel_spec(self,kernel_name):
            return KernelSpec(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],
                              display_name='Python 3',language='python')
    path=ROOT/'notebooks/Double_Ore_14641_PI_verification.ipynb'
    notebook=nbformat.read(path,as_version=4)
    started=datetime.now(timezone.utc).isoformat()
    counts=sum(c.cell_type=='code' for c in notebook.cells)
    for c in notebook.cells:
        if c.cell_type=='code':c.outputs=[];c.execution_count=None
    with tempfile.TemporaryDirectory(prefix='ore14641-kernel-') as clean:
        if kernel_mode=='inprocess':
            from ipykernel.inprocess.manager import InProcessKernelManager
            from ipykernel.inprocess.blocking import BlockingInProcessKernelClient
            class CompatibleClient(BlockingInProcessKernelClient):
                def wait_for_ready(self,timeout=None):
                    # The official in-process client predates this keyword.
                    return super().wait_for_ready()
                def _dispatch_to_kernel(self,msg):
                    # Bind stream output to the request before temporary IO redirection.
                    self.kernel.stdout.set_parent(msg)
                    self.kernel.stderr.set_parent(msg)
                    return super()._dispatch_to_kernel(msg)
            class CompatibleManager(InProcessKernelManager):
                def client(self,**kwds):return CompatibleClient(kernel=self.kernel,session=self.session)
                def shutdown_kernel(self,**kwds):return super().shutdown_kernel()
            km=CompatibleManager()
        else:
            km=KernelManager(kernel_name='python3',kernel_spec_manager=CurrentPython())
        client=NotebookClient(notebook,km=km,timeout=600,kernel_name='python3',
                             resources={'metadata':{'path':clean}},record_timing=True)
        print(f'Executing {counts} code cells in a new IPython kernel ({kernel_mode}); isolated working directory.',flush=True)
        previous_cwd=Path.cwd()
        os.chdir(clean)
        try:
            client.execute()
        except BaseException:
            nbformat.write(notebook,path)
            raise
        finally:
            os.chdir(previous_cwd)
            if kernel_mode=='inprocess' and km.has_kernel:km.shutdown_kernel()
        errors=[output for c in notebook.cells if c.cell_type=='code' for output in c.outputs if output.output_type=='error']
        if errors:raise RuntimeError('Jupyter returned error outputs')
        produced=json.loads((Path(clean)/'notebook_reports/verification.json').read_text())
        independent=json.loads((Path(clean)/'independent_report.json').read_text())
        execution={'method':'nbclient; real ipykernel '+kernel_mode+'; empty temporary working directory',
            'jupyter_kernel_executed':True,'kernel_mode':kernel_mode,
            'kernel_transport':'in-process Jupyter messages' if kernel_mode=='inprocess' else 'local ZMQ',
            'external_kernel_process_tested':kernel_mode=='process','colab_runtime_executed':False,
            'code_cells_executed':counts,'execution_counts':[
                c.execution_count for c in notebook.cells if c.cell_type=='code'],
            'errors':len(errors),'repository_imports_required':False,
            'verification_checks':produced['summary']['total'],
            'independent_checks':independent['passed'],'python':platform.python_version(),
            'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),
            'embedded_source_sha256':{c.metadata['source_file']:c.metadata['source_sha256']
                for c in notebook.cells if 'source_file' in c.metadata}}
        notebook.metadata['execution']=execution
        nbformat.write(notebook,path)
        execution['executed_notebook_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        (ROOT/'reports/notebook_execution.json').write_text(json.dumps(execution,indent=2)+'\n')
        print(f'Executed {counts} cells: {produced["summary"]["total"]} primary and {independent["passed"]} independent checks; zero errors.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--kernel',choices=('process','inprocess'),default='process',
                        help='The official in-process IPython kernel works where local socket binding is unavailable.')
    execute(parser.parse_args().kernel)
