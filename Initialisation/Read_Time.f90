! ============================================================================
! Name        : Read_Time.f90
! Author      : Tristan Salles
! Copyright (C) 2014
! ============================================================================
!> \file Read_Time.f90
!
! Description : Read_Time is used to gather the time information used in badlands.
!
!<
! ============================================================================
module xml_time

    use parallel
    use FoX_sax
    use file_data
    use time_data
    use FoX_common

    implicit none

    public

    ! Time Parameters
    logical, save :: timesection = .false.
    logical, save :: in_Time = .false.
    logical, save :: in_startTime = .false.
    logical, save :: in_endTime = .false.
    logical, save :: in_samplingInterval = .false.
    logical, save :: in_rainInterval = .false.
    logical, save :: in_flowInterval = .false.
    logical, save :: in_displayTime = .false.
    logical, save :: in_outputTime = .false.
    character(len=128), save :: startTime,endTime,outputTime, rainInterval
    character(len=128), save :: displayTime,samplingInterval,flowInterval

contains

    ! ============================================================================
    subroutine startDocument_handler

    end subroutine startDocument_handler
    ! ============================================================================
    subroutine startElement_handler(namespaceURI, localname, name, atts)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name
        type(dictionary_t), intent(in) :: atts

        ! Time Element
        if (name=='Time') in_Time = .true.
        if(in_Time) timesection = .true.
        if(in_Time) call StimeElement_handler(name)

    end subroutine startElement_handler

    !============================================================================
    subroutine endElement_handler(namespaceURI, localname, name)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name

        ! Time Element
        call EtimeElement_handler(name)

    end subroutine endElement_handler
    ! ============================================================================
    subroutine characters_handler(chars)

        character(len=*), intent(in) :: chars

        ! Get Time Parameters
        if (in_Time) call time_characters_handler(chars)

    end subroutine characters_handler
    ! ============================================================================
    subroutine StimeElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='startTime') in_startTime = .true.
        if (name=='endTime') in_endTime = .true.
        if (name=='layerTime') in_displayTime = .true.
        if (name=='outputTime') in_outputTime = .true.
        if (name=='flowInterval') in_flowInterval = .true.
        if (name=='rainInterval') in_rainInterval = .true.
        if (name=='samplingInterval') in_samplingInterval = .true.

    end subroutine StimeElement_handler
    ! ============================================================================
    subroutine EtimeElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='Time') in_Time = .false.
        if (name=='startTime') in_startTime = .false.
        if (name=='endTime') in_endTime = .false.
        if (name=='layerTime') in_displayTime = .false.
        if (name=='outputTime') in_outputTime = .false.
        if (name=='flowInterval') in_flowInterval = .false.
        if (name=='rainInterval') in_rainInterval = .false.
        if (name=='samplingInterval') in_samplingInterval = .false.

    end subroutine EtimeElement_handler
    ! ============================================================================
    subroutine time_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_startTime) then
            startTime = chars
            call rts(startTime,time_start)
        elseif(in_endTime)then
            endTime = chars
            call rts(endTime,time_end)
        elseif(in_outputTime)then
            outputTime = chars
            call rts(outputTime,time_output)
        elseif(in_displayTime)then
            displayTime = chars
            call rts(displayTime,time_layer)
        elseif(in_rainInterval)then
            rainInterval = chars
            call rts(rainInterval,rain_int)
        elseif(in_flowInterval)then
            flowInterval = chars
            call rts(flowInterval,flow_int)
        elseif(in_samplingInterval)then
            samplingInterval = chars
            call rts(samplingInterval,sampling_int)
        endif

    end subroutine time_characters_handler
    ! ============================================================================
    subroutine endDocument_handler

    end subroutine endDocument_handler
    ! ============================================================================
    !> Subroutine xml_reader_time
    !! reads input that defines the time for SP Model experiment (XmL)
    !<
    ! ============================================================================
    subroutine xml_reader_time

        type(xml_t) :: xf
        integer :: ios

        ! Open file
        call open_xml_file(xf, finput , ios)
        if (ios/=0) then
            print*,'---------------------'
            print*, 'Error opening input file for parsing XmL'
            print*,'---------------------'
        endif
        ! Parse file at first and find data
        call parse(xf, &
            startDocument_handler = startDocument_handler, &
            startElement_handler = startElement_handler, &
            endElement_handler = endElement_handler, &
            characters_handler = characters_handler, &
            endDocument_handler = endDocument_handler &
            )
        ! Close file
        call close_xml_t(xf)

    end subroutine xml_reader_time
    ! ============================================================================

end module xml_time
