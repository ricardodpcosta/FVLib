! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: dictionaries
! Modification: February, 2025

! ----------------------------------------------------------------------------
! meshes dictionary
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Reading meshes dictionary")

! check meshes path option
if(fvl_findopt("-m",meshesdictpath) .or. fvl_findopt("--meshes",meshesdictpath)) then
      if(meshesdictpath==fvl_char_empty) then
            call fvl_logerror("Missing argument for option.")
      end if
else
      meshesdictpath          = "setup/meshes.fvd"
end if

! read dictionary
meshesdictfile          = fvl_dict_file(path=meshesdictpath,form="ascii",action="read")
call meshesdictfile     % open()
call meshesdictfile     % readdicts()
call meshesdictfile     % close()

! ----------------------------------------------------------------------------
! model dictionary
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Reading models dictionary")

! check models path option
if(fvl_findopt("-o",modelsdictpath) .or. fvl_findopt("--models",modelsdictpath)) then
      if(modelsdictpath==fvl_char_empty) then
            call fvl_logerror("Missing argument for option.")
      end if
else
      modelsdictpath          = "setup/models.fvd"
end if

! read dictionary
modelsdictfile          = fvl_dict_file(path=modelsdictpath,form="ascii",action="read")
call modelsdictfile     % open()
call modelsdictfile     % readdicts()
call modelsdictfile     % close()

! ----------------------------------------------------------------------------
! schemes dictionary
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Reading schemes dictionary")

! check schemes path option
if(fvl_findopt("-s",schemesdictpath) .or. fvl_findopt("--schemes",schemesdictpath)) then
      if(schemesdictpath==fvl_char_empty) then
            call fvl_logerror("Missing argument for option.")
      end if
else
      schemesdictpath         = "setup/schemes.fvd"
end if

! read dictionary
schemesdictfile         = fvl_dict_file(path=schemesdictpath,form="ascii",action="read")
call schemesdictfile    % open()
call schemesdictfile    % readdicts()
call schemesdictfile    % close()

! ----------------------------------------------------------------------------
! solution dictionary
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Reading solution dictionary")

! check solution path option
if(fvl_findopt("-l",solutiondictpath) .or. fvl_findopt("--solution",solutiondictpath)) then
      if(solutiondictpath==fvl_char_empty) then
            call fvl_logerror("Missing argument for option.")
      end if
else
      solutiondictpath        = "setup/solution.fvd"
end if

! read dictionary
solutiondictfile        = fvl_dict_file(path=solutiondictpath,form="ascii",action="read")
call solutiondictfile   % open()
call solutiondictfile   % readdicts()
call solutiondictfile   % close()

! ----------------------------------------------------------------------------
! control dictionary
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Reading control dictionary")

! check control path option
if(fvl_findopt("-c",controldictpath) .or. fvl_findopt("--control",controldictpath)) then
      if(controldictpath==fvl_char_empty) then
            call fvl_logerror("Missing argument for option -c.")
      end if
else
      controldictpath         = "setup/control.fvd"
end if

! read dictionary
controldictfile         = fvl_dict_file(path=controldictpath,form="ascii",action="read")
call controldictfile    % open()
call controldictfile    % readdicts()
call controldictfile    % close()

! read control parameters
solutionfiledirectory   = controldictfile%getvalue("solution","file_directory",fvl_char_empty)
solutionfileform        = controldictfile%getvalue("solution","file_form",fvl_char_empty)
solutiontimeprecision   = controldictfile%getvalue("solution","time_precision",8)

! ----------------------------------------------------------------------------
! post dictionary
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Reading post dictionary")

! check post path option
if(fvl_findopt("-c",postdictpath) .or. fvl_findopt("--post",postdictpath)) then
      if(postdictpath==fvl_char_empty) then
            call fvl_logerror("Missing argument for option -c.")
      end if
else
      postdictpath            = "setup/post.fvd"
end if

! read dictionary
postdictfile            = fvl_dict_file(path=postdictpath,form="ascii",action="read")
call postdictfile       % open()
call postdictfile       % readdicts()
call postdictfile       % close()

! end of file
