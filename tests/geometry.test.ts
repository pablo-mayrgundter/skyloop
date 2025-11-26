import { Vector3 } from 'three'
import { EARTH_RADIUS_KM, LOOP_RADIUS_KM } from '../src/config/simulation'
import { greatCirclePoint, interpolateGreatCircle } from '../src/config/simulation'

describe('geometry utilities', () => {
  describe('distance calculations', () => {
    it('calculates Euclidean distance between two points', () => {
      const p1 = greatCirclePoint(0, 0, EARTH_RADIUS_KM)
      const p2 = greatCirclePoint(0, 90, EARTH_RADIUS_KM)
      
      // For 90° separation on a sphere, chord distance = √2 * radius
      // Arc distance = π/2 * radius, but Euclidean distance is chord
      const expectedChordDistance = Math.sqrt(2) * EARTH_RADIUS_KM
      const actualDistance = p1.distanceTo(p2)
      
      expect(Math.abs(actualDistance - expectedChordDistance) / expectedChordDistance).toBeLessThan(0.01)
    })
    
    it('calculates great circle arc distance', () => {
      const p1 = greatCirclePoint(0, 0, EARTH_RADIUS_KM)
      const p2 = greatCirclePoint(0, 90, EARTH_RADIUS_KM)
      
      // Calculate arc distance using dot product
      const dot = p1.clone().normalize().dot(p2.clone().normalize())
      const angle = Math.acos(Math.max(-1, Math.min(1, dot)))
      const arcDistance = angle * EARTH_RADIUS_KM
      
      // 90° = π/2 radians
      const expectedArcDistance = (Math.PI / 2) * EARTH_RADIUS_KM
      expect(Math.abs(arcDistance - expectedArcDistance) / expectedArcDistance).toBeLessThan(0.01)
    })

    it('calculates distance between LA and SF', () => {
      const la = greatCirclePoint(34.0522, -118.2437, EARTH_RADIUS_KM)
      const sf = greatCirclePoint(37.7749, -122.4194, EARTH_RADIUS_KM)
      
      const distance = la.distanceTo(sf)
      
      // LA to SF is approximately 560km
      expect(distance).toBeGreaterThan(500) // 500km
      expect(distance).toBeLessThan(600) // 600km
    })
  })

  describe('vector operations', () => {
    it('normalizes vectors correctly', () => {
      const point = greatCirclePoint(45, 90, EARTH_RADIUS_KM)
      const normalized = point.clone().normalize()
      
      expect(normalized.length()).toBeCloseTo(1, 5)
    })

    it('calculates dot product for angle', () => {
      const p1 = greatCirclePoint(0, 0, EARTH_RADIUS_KM)
      const p2 = greatCirclePoint(0, 90, EARTH_RADIUS_KM)
      const p3 = greatCirclePoint(0, 180, EARTH_RADIUS_KM)
      
      const dot1 = p1.clone().normalize().dot(p2.clone().normalize())
      const dot2 = p1.clone().normalize().dot(p3.clone().normalize())
      
      // p1 and p2 are 90° apart, dot should be ~0
      expect(Math.abs(dot1)).toBeLessThan(0.1)
      
      // p1 and p3 are 180° apart, dot should be ~-1
      expect(dot2).toBeCloseTo(-1, 1)
    })

    it('calculates cross product', () => {
      const p1 = greatCirclePoint(0, 0, EARTH_RADIUS_KM)
      const p2 = greatCirclePoint(0, 90, EARTH_RADIUS_KM)
      
      const cross = new Vector3().crossVectors(
        p1.clone().normalize(),
        p2.clone().normalize()
      )
      
      // Cross product should be perpendicular to both
      expect(cross.length()).toBeGreaterThan(0)
    })
  })

  describe('interpolation properties', () => {
    it('maintains constant speed along arc', () => {
      const steps = 10
      const distances: number[] = []
      
      for (let i = 0; i < steps; i++) {
        const t1 = i / steps
        const t2 = (i + 1) / steps
        const p1 = interpolateGreatCircle(
          { lat: 0, lon: 0 },
          { lat: 0, lon: 90 },
          t1,
          EARTH_RADIUS_KM
        )
        const p2 = interpolateGreatCircle(
          { lat: 0, lon: 0 },
          { lat: 0, lon: 90 },
          t2,
          EARTH_RADIUS_KM
        )
        distances.push(p1.distanceTo(p2))
      }
      
      // Distances should be roughly equal (allowing for floating point)
      const avgDistance = distances.reduce((a, b) => a + b, 0) / distances.length
      for (const dist of distances) {
        expect(Math.abs(dist - avgDistance)).toBeLessThan(avgDistance * 0.2) // 20% tolerance
      }
    })

    it('handles interpolation at loop radius', () => {
      const point = interpolateGreatCircle(
        { lat: 34.0522, lon: -118.2437 },
        { lat: 37.7749, lon: -122.4194 },
        0.5,
        LOOP_RADIUS_KM
      )
      
      expect(Math.abs(point.length() - LOOP_RADIUS_KM)).toBeLessThan(1000)
    })
  })

  describe('coordinate system', () => {
    it('uses right-handed coordinate system', () => {
      const equator = greatCirclePoint(0, 0, EARTH_RADIUS_KM)
      const equator90 = greatCirclePoint(0, 90, EARTH_RADIUS_KM)
      
      // Cross product should give correct orientation
      const cross = new Vector3().crossVectors(
        equator.clone().normalize(),
        equator90.clone().normalize()
      )
      
      // Should point in positive Y direction (north pole)
      expect(cross.y).toBeGreaterThan(0)
    })

    it('handles all quadrants correctly', () => {
      const points = [
        greatCirclePoint(45, 45, EARTH_RADIUS_KM),   // NE
        greatCirclePoint(45, -45, EARTH_RADIUS_KM),  // NW
        greatCirclePoint(-45, 45, EARTH_RADIUS_KM),  // SE
        greatCirclePoint(-45, -45, EARTH_RADIUS_KM)  // SW
      ]
      
      for (const point of points) {
        expect(point.length()).toBeCloseTo(EARTH_RADIUS_KM, -3)
      }
    })
  })

  describe('edge cases and boundaries', () => {
    it('handles extreme latitudes', () => {
      const nearPole = greatCirclePoint(89.9, 0, EARTH_RADIUS_KM)
      expect(nearPole.length()).toBeCloseTo(EARTH_RADIUS_KM, -3)
    })

    it('handles extreme longitudes', () => {
      const point1 = greatCirclePoint(0, 179.9, EARTH_RADIUS_KM)
      const point2 = greatCirclePoint(0, -179.9, EARTH_RADIUS_KM)
      expect(point1.length()).toBeCloseTo(EARTH_RADIUS_KM, -3)
      expect(point2.length()).toBeCloseTo(EARTH_RADIUS_KM, -3)
    })

    it('handles zero radius', () => {
      const point = greatCirclePoint(45, 90, 0)
      expect(point.length()).toBe(0)
    })

    it('handles very large radius', () => {
      const largeRadius = 10000000 // 10,000km
      const point = greatCirclePoint(45, 90, largeRadius)
      expect(point.length()).toBeCloseTo(largeRadius, -2)
    })
  })
})

