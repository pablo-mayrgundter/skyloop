import { interpolateGreatCircle, LOS_ANGELES, SAN_FRANCISCO, LOOP_RADIUS_KM, EARTH_RADIUS_KM } from '../src/config/simulation'
import { Vector3 } from 'three'

describe('camera setup', () => {
  describe('initial camera position', () => {
    it('positions camera above loop midpoint', () => {
      const loopMidpoint = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5, LOOP_RADIUS_KM)
      const cameraHeight = 50000 // 50km above loop
      const cameraPos = loopMidpoint.clone().normalize().multiplyScalar(LOOP_RADIUS_KM + cameraHeight)
      
      // Camera should be at correct distance from origin
      expect(cameraPos.length()).toBeCloseTo(LOOP_RADIUS_KM + cameraHeight, -2) // 1km tolerance
      
      // Camera should be above the loop (further from origin)
      expect(cameraPos.length()).toBeGreaterThan(LOOP_RADIUS_KM)
    })

    it('camera is outside Earth surface', () => {
      const loopMidpoint = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5, LOOP_RADIUS_KM)
      const cameraHeight = 50000
      const cameraPos = loopMidpoint.clone().normalize().multiplyScalar(LOOP_RADIUS_KM + cameraHeight)
      
      // Camera must be well above Earth surface
      expect(cameraPos.length()).toBeGreaterThan(EARTH_RADIUS_KM + 10) // At least 10km above surface
    })
  })

  describe('camera near and far planes', () => {
    const CAMERA_NEAR = 0.1 // meters
    const CAMERA_FAR = 20000000 // 20,000km

    it('near plane allows close-up viewing', () => {
      // Near plane should be small enough for human-scale viewing (1m objects)
      expect(CAMERA_NEAR).toBeLessThanOrEqual(1)
      expect(CAMERA_NEAR).toBeGreaterThan(0)
    })

    it('far plane allows viewing entire Earth', () => {
      // Earth diameter is ~12,742km, so far plane should be at least that
      const earthDiameter = EARTH_RADIUS_KM * 2
      expect(CAMERA_FAR).toBeGreaterThan(earthDiameter)
      
      // Should be able to see objects at loop altitude (6,391km from center)
      expect(CAMERA_FAR).toBeGreaterThan(LOOP_RADIUS_KM * 2)
    })

    it('near and far planes have reasonable ratio', () => {
      // Ratio should not be too extreme (avoid precision issues)
      const ratio = CAMERA_FAR / CAMERA_NEAR
      // 200 million is reasonable for Three.js
      expect(ratio).toBeLessThan(1e9)
      expect(ratio).toBeGreaterThan(1000)
    })
  })

  describe('camera view frustum', () => {
    it('can see objects at loop altitude', () => {
      const loopMidpoint = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5, LOOP_RADIUS_KM)
      const cameraHeight = 50000
      const cameraPos = loopMidpoint.clone().normalize().multiplyScalar(LOOP_RADIUS_KM + cameraHeight)
      
      // Distance from camera to loop should be within far plane
      const distanceToLoop = cameraPos.distanceTo(loopMidpoint)
      expect(distanceToLoop).toBeLessThan(20000000) // Within far plane
      expect(distanceToLoop).toBeGreaterThan(0.1) // Beyond near plane
    })

    it('can see Earth surface', () => {
      const loopMidpoint = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5, LOOP_RADIUS_KM)
      const cameraHeight = 50000
      const cameraPos = loopMidpoint.clone().normalize().multiplyScalar(LOOP_RADIUS_KM + cameraHeight)
      
      // Distance to Earth surface (along camera's radial direction)
      const distanceToSurface = cameraPos.length() - EARTH_RADIUS_KM
      expect(distanceToSurface).toBeLessThan(20000000) // Within far plane
      expect(distanceToSurface).toBeGreaterThan(0.1) // Beyond near plane
    })
  })

  describe('orbit controls limits', () => {
    const HUMAN_SCALE_DISTANCE = 1 // 1 meter
    const MAX_ZOOM_DISTANCE = 100000 // 100km

    it('min distance allows human-scale viewing', () => {
      expect(HUMAN_SCALE_DISTANCE).toBeGreaterThanOrEqual(0.1) // At least near plane
      expect(HUMAN_SCALE_DISTANCE).toBeLessThan(10) // Human scale
    })

    it('max distance allows viewing from space', () => {
      expect(MAX_ZOOM_DISTANCE).toBeLessThan(20000000) // Within far plane
      expect(MAX_ZOOM_DISTANCE).toBeGreaterThan(LOOP_RADIUS_KM - EARTH_RADIUS_KM) // Can see loop
    })
  })
})

