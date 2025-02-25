# Setting up FoundHaplo MySQL Database

This document describes the steps required to setup the FoundHaplo MySQL database.


Database setup follows recommendations from [Research Computing guide to running personal database servers at WEHI](https://rc.wehi.edu.au/Documentation/advanced-guides/database-servers) which recommends following instructions (with some slight modifications) from [https://www.hpc.iastate.edu/guides/containers/mysql-server](https://www.hpc.iastate.edu/guides/containers/mysql-server)

## How to initiate the database

1. Create a folder for the FoundHaplo database

```bash
mkdir -p /mypath/FoundHaplo_database
FoundHaplo_database_DIR=/mypath/FoundHaplo_database
cd FoundHaplo_database_DIR
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
The MySQL server is now running on the server

9. To connect to the server we use the installed `mysql` client with the additional information about where to find the local socket (directory bound to container)

```bash
mysql -S $FoundHaplo_database_DIR/mysql/run/mysqld/mysqld.sock
```
You are now running MySQL server using singularity!

10. Create the FoundHaplo database (FoundHaploDB) and tables in MySQL using schema as below.

```
source /mypath/FoundHaplo/FoundHaplo_database/FoundHaplo_database_create_schema.sql;
exit;
```
11. Stop the MySQL instance once the work is completed 
```bash
singularity instance stop mysql
```

## How to acess the database

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
2. To create a remote user to access the database by multiple nodes if needed. Note down the generated password 
```bash
singularity exec instance://mysql create_remote_admin_user.sh 
```

3. Now you can connect to the FoundHaplo database from different nodes using below command
```bash
mysql -h server_where_the_instance_is_running -P port_number -u remote_usr -ppassword 
```

4. Stop the MySQL instance once the work is completed 
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
export PATH=$PATH:/wehisan/bioinf/lab_bahlo/projects/methods_dev/Haplotype_Dee/Software/mariadb-connector-c-3.1.11-centos7-amd64/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/wehisan/bioinf/lab_bahlo/projects/methods_dev/Haplotype_Dee/Software/mariadb-connector-c-3.1.11-centos7-amd64/lib/mariadb
```

3. install RMariaDB in R
```R
install.packages("RMariaDB") 
```

4. Run /mypath/FoundHaplo/FoundHaplo_database/Create_SQL_script_to_import.R to import each created disease haplotype into the database.
```R
Create_SQL_script_to_import(disease_hap_file,save_SQL_file,port,host,password,dbname,unix.socket,family_id,individual_id,father_id,mother_id,sex,ethnicity,ethnicity_superpopulation,ethnicity_method,sample_id,data_type,external_lab_id,impute_method,impute_panel,import_date,mutation_id,disease,disease_id,omim_id,gene=inheritance_model,chr,start_hg19,end_hg19,start_hg38,end_hg38,start_cm,end_cm,genotype,validated,validation_method,validation_note,format_vcf)
```

5. Connect to the database and source the resulting sql script into the FoundHaploDB to import disease haplotypes

Go back to the [documentaton](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).


