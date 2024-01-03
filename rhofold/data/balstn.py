# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import subprocess
from typing import Any, Mapping
from rhofold.utils import tmpdir, timing
import shutil, os

class BLASTN:
  """Python wrapper of the BLAST binary."""

  def __init__(self,
               *,
               binary_dpath: str,
               databases,
               n_cpu: int = 4,
               ):
    """Initializes the Python BLAST wrapper.

    Args:
      binary_dpath: The path to the BLAST executable.
      databases: A sequence of BLAST database paths.
    """

    self.binary_path = binary_dpath
    self.databases = databases
    self.n_cpu = n_cpu

  def query(self, input_fasta_path, output_msa_path, logger) -> Mapping[str, Any]:
      """Queries the database using BLAST."""


      with tmpdir(base_dir=os.path.dirname(output_msa_path)) as query_tmp_dir:

          blast_aln = f'{query_tmp_dir}/blast.aln'
          a3m_path = os.path.join(query_tmp_dir, f'blast.a3m')

          blast_outputs = []

          for database in self.databases:
              db_name = os.path.basename(database)
              blast = f'{query_tmp_dir}/blast_{db_name}.db'

              cmd = [f'{self.binary_path}/blastn',
                    '-db', database,
                    '-query',input_fasta_path,
                    '-out', blast,
                    '-evalue', '0.001',
                    '-num_descriptions','1',
                    '-num_threads', f'{self.n_cpu}',
                    '-line_length', '1000',
                    '-num_alignments','50000',
                    '-strand','plus','-task','blastn'
                     ]

              process = subprocess.Popen(
                  cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

              with timing(f'BLASTN scearch @ {db_name}', logger=logger):
                stdout, stderr = process.communicate()
                retcode = process.wait()

              if retcode:
                logger.error('    BLAST failed. BLAST stderr begin:')
                for error_line in stderr.decode('utf-8').splitlines():
                  if error_line.strip():
                    logger.error(error_line.strip())
                logger.error('    BLAST stderr end')
                raise RuntimeError('BLAST failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                    stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

              blast_outputs.append(blast)

          blast = f'{query_tmp_dir}/blast.db'
          cmd = f"cat {' '.join(blast_outputs)} > {blast} "
          os.system(cmd)

          cmd = f"perl {self.binary_path}/parse_blastn_local.pl {blast} {input_fasta_path} {blast_aln} > {query_tmp_dir}/log.txt \n" \
                f"sed -i 's/\s.*$//' {blast_aln} \n" \
                f"sed -i 's/$seq_id/$seq_id E=0.0/g' {blast_aln} \n" \
                f"perl {self.binary_path}/reformat.pl fas a3m {blast_aln} {a3m_path} > {query_tmp_dir}/log.txt"
          os.system(cmd)

          shutil.copy(a3m_path, output_msa_path)


