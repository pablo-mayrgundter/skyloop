import { Vector3 } from 'three'
import {
  EARTH_RADIUS_KM,
  LOOP_RADIUS_KM,
  LOS_ANGELES,
  SAN_FRANCISCO,
  greatCirclePoint,
  interpolateGreatCircle
} from '../src/config/simulation'

describe('greatCirclePoint', () => {
  it('returns point at correct radius', () => {
    const la = greatCirclePoint(LOS_ANGELES.lat, LOS_ANGELES.lon, EARTH_RADIUS_KM)
    expect(Math.abs(la.length() - EARTH_RADIUS_KM)).toBeLessThan(1000) // 1km tolerance for floating point
  })

  it('handles different radii', () => {
    const customRadius = 1000000 // 1000km
    const point = greatCirclePoint(0, 0, customRadius)
    expect(Math.abs(point.length() - customRadius)).toBeLessThan(100)
  })

  it('returns correct coordinates for North Pole', () => {
    const northPole = greatCirclePoint(90, 0, EARTH_RADIUS_KM)
    expect(Math.abs(northPole.x)).toBeLessThan(1000) // Within 1km of 0
    expect(Math.abs(northPole.z)).toBeLessThan(1000)
    expect(Math.abs(northPole.y - EARTH_RADIUS_KM)).toBeLessThan(1000)
  })

  it('returns correct coordinates for South Pole', () => {
    const southPole = greatCirclePoint(-90, 0, EARTH_RADIUS_KM)
    expect(Math.abs(southPole.x)).toBeLessThan(1000)
    expect(Math.abs(southPole.z)).toBeLessThan(1000)
    expect(Math.abs(southPole.y + EARTH_RADIUS_KM)).toBeLessThan(1000)
  })

  it('returns correct coordinates for equator at 0° longitude', () => {
    const equator = greatCirclePoint(0, 0, EARTH_RADIUS_KM)
    expect(Math.abs(equator.y)).toBeLessThan(1000) // Should be near 0
    // x should be approximately -EARTH_RADIUS_KM (negative due to coordinate system)
    expect(Math.abs(Math.abs(equator.x) - EARTH_RADIUS_KM)).toBeLessThan(1000)
    expect(Math.abs(equator.z)).toBeLessThan(1000) // Should be near 0
    expect(equator.length()).toBeCloseTo(EARTH_RADIUS_KM, -3) // Radius should be correct
  })

  it('handles longitude wrapping', () => {
    const point1 = greatCirclePoint(0, 180, EARTH_RADIUS_KM)
    const point2 = greatCirclePoint(0, -180, EARTH_RADIUS_KM)
    expect(point1.distanceTo(point2)).toBeLessThan(1000) // Should be same point
  })

  it('is consistent for same coordinates', () => {
    const p1 = greatCirclePoint(45, 90, EARTH_RADIUS_KM)
    const p2 = greatCirclePoint(45, 90, EARTH_RADIUS_KM)
    expect(p1.distanceTo(p2)).toBeLessThan(0.001)
  })
})

describe('interpolateGreatCircle', () => {
  it('interpolates between endpoints', () => {
    const mid = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5)
    expect(mid.length()).toBeCloseTo(EARTH_RADIUS_KM, -3) // 1km tolerance
    const start = greatCirclePoint(LOS_ANGELES.lat, LOS_ANGELES.lon)
    const end = greatCirclePoint(SAN_FRANCISCO.lat, SAN_FRANCISCO.lon)
    expect(mid.distanceTo(start)).toBeGreaterThan(0)
    expect(mid.distanceTo(end)).toBeGreaterThan(0)
  })

  it('returns start point at t=0', () => {
    const start = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0)
    const expectedStart = greatCirclePoint(LOS_ANGELES.lat, LOS_ANGELES.lon)
    expect(start.distanceTo(expectedStart)).toBeLessThan(1000) // 1km tolerance
  })

  it('returns end point at t=1', () => {
    const end = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 1)
    const expectedEnd = greatCirclePoint(SAN_FRANCISCO.lat, SAN_FRANCISCO.lon)
    expect(end.distanceTo(expectedEnd)).toBeLessThan(1000) // 1km tolerance
  })

  it('maintains radius for all interpolation values', () => {
    for (let t = 0; t <= 1; t += 0.1) {
      const point = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t)
      expect(Math.abs(point.length() - EARTH_RADIUS_KM)).toBeLessThan(1000)
    }
  })

  it('handles custom radius', () => {
    const customRadius = LOOP_RADIUS_KM
    const point = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5, customRadius)
    expect(Math.abs(point.length() - customRadius)).toBeLessThan(1000)
  })

  it('handles same start and end point', () => {
    const point = interpolateGreatCircle(LOS_ANGELES, LOS_ANGELES, 0.5)
    const expected = greatCirclePoint(LOS_ANGELES.lat, LOS_ANGELES.lon)
    expect(point.distanceTo(expected)).toBeLessThan(1000)
  })

  it('interpolates smoothly (monotonic distance)', () => {
    const start = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0)
    const end = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 1)
    const distances: number[] = []
    
    for (let t = 0; t <= 1; t += 0.1) {
      const point = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t)
      distances.push(start.distanceTo(point))
    }
    
    // Distance should generally increase (allowing for some floating point noise)
    for (let i = 1; i < distances.length; i++) {
      expect(distances[i]).toBeGreaterThanOrEqual(distances[i - 1] - 1000) // Allow 1km tolerance
    }
  })

  it('handles edge cases: t < 0', () => {
    const point = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, -0.1)
    expect(point.length()).toBeCloseTo(EARTH_RADIUS_KM, -3)
  })

  it('handles edge cases: t > 1', () => {
    const point = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 1.1)
    expect(point.length()).toBeCloseTo(EARTH_RADIUS_KM, -3)
  })

  it('produces points on great circle arc', () => {
    const start = greatCirclePoint(LOS_ANGELES.lat, LOS_ANGELES.lon)
    const end = greatCirclePoint(SAN_FRANCISCO.lat, SAN_FRANCISCO.lon)
    const mid = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5)
    
    // Midpoint should be approximately equidistant from start and end
    const distToStart = mid.distanceTo(start)
    const distToEnd = mid.distanceTo(end)
    const totalDist = start.distanceTo(end)
    
    // Midpoint distance should be roughly half the total distance
    expect(Math.abs(distToStart + distToEnd - totalDist)).toBeLessThan(totalDist * 0.1) // 10% tolerance
  })
})
