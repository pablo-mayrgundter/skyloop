export type SimulationMode = 'firstPerson' | 'birdseye'

export type SkyloopModuleId =
  | 'main'
  | 'sections'
  | 'tethers'
  | 'tube'
  | 'gliders'
  | 'atmosphere'
  | 'weather'
  | 'dayNight'
  | 'bom'
  | 'economics'

export type OverlayFlags = {
  showAtmosphere: boolean
  showWind: boolean
  showTethers: boolean
  showSections: boolean
  showGliders: boolean
  showDayNight: boolean
}

export type AtmosphericLayer = {
  name: string
  altitudeKm: number
  densityKgPerM3: number
}

export type WindProfile = {
  altitudeKm: number
  speedMps: number
  directionDeg: number
}

export type CostBreakdown = {
  capexUsd: number
  opexUsdPerYear: number
  notes?: string
}

export type GreatCirclePoint = {
  lat: number
  lon: number
}

export type SimulationMetrics = {
  velocityMps: number
  altitudeKm: number
  tetherTensionKn: number
  powerMw: number
}
