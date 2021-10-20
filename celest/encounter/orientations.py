

def generate_pointing_profiles(self, groundPos: Any, encInd: np.array,
                                   maneuverTime: int) -> Rotation:
        """Generate satellite rotations for ground tracking.

        This function is intended to take in a single ground location along
        with the windows at which the spacecraft makes imaging passes over the
        location. This method uses the Odyssey Pointing Profile
        determination system created by Mingde Yin.

        NOTE: This should be generalized in the future to many sites.

        Parameters
        ----------
        groundPos : GroundPosition
            Ground location and encounter information.
        encInd: np.array
            Array of arrays of indices that correspond to times and positions
            where the spacecraft is in an imaging encounter window with the
            given ground location.
        maneuverTime: float
            Number of array indices to pad on either side of
            an encounter window to use for maneuvering time.
            TODO: Turn into a standard time unit since setting a fixed number
            of array indices is limited.

        Notes
        -----
        The strategy for pointing profile generation is as follows:
        1. The default orientation is to have the spacecraft camera pointing
        towards its zenith.
        2. When the spacecraft is imaging, orient the satellite such that the
        camera is facing the target.
        3. Generate an initial set of pointing profiles assuming the above.
        4. Interpolate the rotations between the normal and target-acquired
           states to smooth out transitions.
        """

        ground_GEO = [groundPos.coor[0], groundPos.coor[1], groundPos.radius]
        ground_GEO = np.repeat(np.array([ground_GEO]), self.position.length, axis=0)
        target_site = Coordinate(ground_GEO, "GEO", self.times)

        # Stage 1: preliminary rotations.
        # Get difference vector between spacecraft and target site.
        SC_to_site = target_site.ECI() - self.position.ECI()

        # Point toward the zenith except when over the imaging site.
        pointing_directions = self.position.ECI()
        pointing_directions[encInd, :] = SC_to_site[encInd, :]

        # Preliminary rotation set.
        # Temporarily represent as quaternion for interpolation.
        rotations = sat_rotation(pointing_directions).as_quat()

        # Set flight modes. By default, point normal.
        flight_ind = 0 * np.ones(SC_to_site.shape[0])

        # Point to target during encounters.
        flight_ind[encInd] = 2

        # Stage 2: interpolation.
        # Generate sets of encounters which are clustered together.
        # This takes individual encInd into clusters which we can use later to
        # figure out when to start interpolation.
        split_ind = np.where(np.diff(encInd) > 1)[0]+1
        encounter_segments = np.split(encInd, split_ind)

        for encInd in encounter_segments:

            start_step = encInd[0] - maneuverTime
            end_step = encInd[-1] + maneuverTime

            # Get starting and ending quaternions.
            start_rotation = rotations[start_step]
            end_rotation = rotations[end_step]

            slerp_1 = Slerp([start_step, encInd[0]], Rotation.from_quat([start_rotation, rotations[encInd[0]]]))
            interp_rotations_1 = slerp_1(np.arange(start_step, encInd[0]))

            rotations[start_step:encInd[0]] = interp_rotations_1.as_quat()
            flight_ind[start_step:encInd[0]] = 1

            slerp_2 = Slerp([encInd[-1], end_step], Rotation.from_quat([rotations[encInd[-1]], end_rotation]))
            interp_rotations_2 = slerp_2(np.arange(encInd[-1], end_step))

            rotations[encInd[-1]:end_step] = interp_rotations_2.as_quat()
            flight_ind[encInd[-1]:end_step] = 3

        return Rotation.from_quat(rotations)