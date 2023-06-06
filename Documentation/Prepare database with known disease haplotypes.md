# Setting up FoundHaplo MySQL Database

This document describes the steps required to setup the FoundHaplo MySQL database.


Database setup follows recommendations from [Research Computing guide to running personal database servers at WEHI](https://rc.wehi.edu.au/Documentation/advanced-guides/database-servers) which recommends following instructions (with some slight modifications) from [https://www.hpc.iastate.edu/guides/containers/mysql-server](https://www.hpc.iastate.edu/guides/containers/mysql-server)

## How to initiate the database

1. Create a folder for the FoundHaplo database
```bash
mkdir -p /mypath/FoundHaplo_database
FoundHaplo_database_DIR=/mypath/FoundHaplo_database
```
Or create a folder "FoundHaplo_database" for the database inside FoundHaplo_DIR/FoundHaplo_database directory
```bash
mkdir -p $FoundHaplo_DIR/FoundHaplo_database/FoundHaplo_database
FoundHaplo_database_DIR=$FoundHaplo_DIR/FoundHaplo_database/FoundHaplo_database
cd $FoundHaplo_database_DIR
```
2. Setup singularity

Load the singularity module
```bash
module load singularity
```

3. Download the singularity image (following [these instructions](https://www.hpc.iastate.edu/guides/containers/mysql-server) as recommended by the [Research Computing guide](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Database-servers.aspx))

```bash
singularity pull --name mysql.simg shub://ISU-HPC/mysql
mkdir setup
```

4. Download config scripts to update MySQL root password

```bash
curl https://raw.githubusercontent.com/ISU-HPC/mysql/master/my.cnf > setup/my.cnf
curl https://raw.githubusercontent.com/ISU-HPC/mysql/master/mysqlrootpw > setup/mysqlrootpw
```

5. Manually edit `setup/my.cnf`. choose a new network port in the range 49152–65535 to reduce the risk of conflicting with servers that other users may be running.

setup/my.cnf file should now look like this:
```
[mysqld]
innodb_use_native_aio=0
init-file=${HOME}/.mysqlrootpw
port=port in the range 49152–65535
bind-address = 0.0.0.0

[client]
user=root
password='my-secret-pw'
```

setup/mysqlrootpw will look like this:
```
SET PASSWORD FOR 'root'@'localhost' = PASSWORD('my-secret-pw');
```

6. Copy `setup/my.cnf` and `setup/mysqlrootpw` to `~/.my.cnf` and `~/.mysqlrootpw`:

```bash
cp setup/my.cnf ${HOME}/.my.cnf
cp setup/mysqlrootpw ${HOME}/.mysqlrootpw
```

7. Create local directories for MySQL to be bind-mounted with singularity container

```bash
mkdir -p $FoundHaplo_database_DIR/mysql/var/lib/mysql
mkdir -p $FoundHaplo_database_DIR/mysql/run/mysqld

chmod g+w -R $FoundHaplo_database_DIR/mysql
```

8. Start MySQL server

We can now start the singularity instance for the MySQL server

```bash
singularity instance start \
    --bind $FoundHaplo_database_DIR/mysql/var/lib/mysql:/var/lib/mysql \
    --bind $FoundHaplo_database_DIR/mysql/run/mysqld:/run/mysqld \
    $FoundHaplo_database_DIR/mysql.simg mysql
```

```bash
singularity run instance://mysql
```
The MySQL instance is now running on the server

## How to acess the MySQL instance

1. Run singularity instance if its already stopped.
```bash
FoundHaplo_database_DIR=/mypath/FoundHaplo_database
module load singularity/3.6.2
singularity instance stop mysql
cd $FoundHaplo_database_DIR

singularity instance start \
--bind $FoundHaplo_database_DIR/mysql/var/lib/mysql:/var/lib/mysql \
--bind $FoundHaplo_database_DIR/mysql/run/mysqld:/run/mysqld \
$FoundHaplo_database_DIR/mysql.simg mysql

singularity run instance://mysql
```
2. To create a remote user to access the MySQL instance by multiple nodes if needed. Note down the generated password 
```bash
singularity exec instance://mysql create_remote_admin_user.sh 
```

3. Now you can connect to the the MySQL instance from different nodes using below command
```bash
mysql -h server_where_the_instance_is_running -P port_number -u remote_usr -ppassword 
```

You are now connected to MySQL server using singularity!

4. If you have not already created the FoundHaplo database, create the FoundHaplo database (FoundHaploDB) and tables in MySQL using schema as below.

```MySQL
source /mypath/FoundHaplo/FoundHaplo_database/FoundHaplo_database_create_schema.sql;
exit;
```

5. Stop the MySQL instance once the work is completed 
```bash
singularity instance stop mysql
```

## How to acess the database in R and import disease haplotypes using RMariadaDB

1. Download the relevant connector using [link](https://mariadatabase.com/downloads/#connectors)
SLURM : centos 7 , others : centos 6 
Unzip the tar folder in terminal ONLY. Do not unzip by right clicking

2. Add the downloaded files to your default path everytime the database is needed to be accessed from R

Add the mariadatabase config to the $PATH  

Add the directory containing libmariadatabase.so.3 to your LD_LIBRARY_PATH. 

```bash
export PATH=$PATH:/mypath/mariadb-connector-c-3.1.11-centos7-amd64/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mypath/mariadb-connector-c-3.1.11-centos7-amd64/lib/mariadb
```
Create a directory named scripts to save .sql scripts to import disease haplotypes created in R
```bash
mkdir $FoundHaplo_database_DIR/scripts
```

3. install RMariaDB in R
```R
install.packages("RMariaDB") 
```

4. Run [/mypath/FoundHaplo/FoundHaplo_database/Create_SQL_script_to_import.R](https://github.com/bahlolab/FoundHaplo/blob/main/FoundHaplo_database/Create_SQL_script_to_import.R) to import each created disease haplotype into the database. Create_SQL_script_to_import.R script can only import one disease haplotype at a time.
```R
Create_SQL_script_to_import(disease_hap_FILE=/mypath/FoundHaplo/input_files/input_vcf_data/disease_haplotypes/HG00313_1_HG00313_1,HG00313_2_HG00313_2.vcf.gz,save_SQL_FILE=/mypath/FoundHaplo/FoundHaplo_database/FoundHaplo_database/scripts/example.sql,db_port=port_number,db_host=server_where_the_instance_is_running,db_password=pwd,db_name=FoundHaploDB,db_unix_socket=/mypath/FoundHaplo/FoundHaplo_database/FoundHaplo_database/mysql/run/mysqld/mysqld.sock,family_id=1,individual_id=1,father_id=0,mother_id=0,sex=1,sex_method="UNK",ancestral_population="Australian",ancestral_superpopulation="EUR",ancestry_method="reported",sample_id=1,data_type="SNP chip",external_lab_id="external_lab_id_example",external_source="external_source_example",phasing_method="trio phasing",impute_method="MIS",impute_phasing_panel="1000 Genomes EUR hg19",import_date="2023-05-01",DCV_id=1,disease_name="FAME1",omim_id=601068,gene="SAMD12",genomic_region="intronic",inheritance_model="AD",chromosome="chr8",start_position_hg19=19379052,end_position_hg19=-99999,start_position_hg38=-99999,end_position_hg38=-99999,start_position_cM=-99999,end_position_cM=-99999,genotype=1,validated=1,validation_method="RP-PCR,validation_note="UNK")
```
All the parameters that user has to specify are described below

* disease_hap_FILE = Disease haplotype VCF file created using guidelines [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20known%20disease%20haplotypes.md)
* save_SQL_FILE = Path to save the .sql command files 
* db_port= Port number of the database
* db_host= Database host
* db_password= Randomly generated user password to access the database from R
* db_name= Name of the database
* db_unix_socket= Path to .sock of the database i.e $FoundHaplo_database_DIR/mysql/run/mysqld/mysqld.sock

Parameters starting from family_id must be specified based on the databse schema as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/FoundHaplo_database/FoundHaplo_database_info.docx)

5. Connect to the database and source the resulting sql script of the Create_SQL_script_to_import.R saved in save_SQL_FILE into the FoundHaploDB to import disease haplotypes

```MySQL
source /mypath_to_save_SQL_FILE;
exit;
```
Repeat steps 4 and 5 to import all the disease haplotypes into the database

6. Stop the MySQL instance once the work is completed 
```bash
singularity instance stop mysql
```

Note : You can keep connecting to the MySQL instance using "mysql -h server_where_the_instance_is_running -P port_number -u remote_usr -ppassword" as long as the the singularity instance is not terminated. Singularity instance must be terminated after working/accessing the database to help prevent database getting crashed.

Go back to the [documentaton](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).


